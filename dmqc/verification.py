"""Independent, read-only checks for core Argo single-cycle D-profile files.

Coverage: profile format 3.1, User's Manual v3.3a section 2.2 and reference
tables, CTD QC Manual v3.3 sections 3.6 and 4.7. This is a local subset of
GDAC checks, not the official OneArgo format checker or scientific QC.
"""
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset


SOURCES = {
    'user_manual': 'https://cdn.ioos.noaa.gov/media/2020/03/argo_user_manual_v3.3a.pdf',
    'ctd_qc_manual': 'https://cdn.ioos.noaa.gov/media/2020/03/Argo-QC-for-CTD-and-Trajectory-Data.pdf',
    'official_checker': 'https://github.com/OneArgo/ArgoFormatChecker',
}
SCOPE = ('Local core Argo profile-format 3.1 checks; not full GDAC certification. '
         'No climatology, sensor-drift, land-mask, vocabulary-service, BGC or trajectory checks.')
FILE_NAME = re.compile(r'D(?P<float>\d{5,7})_(?P<cycle>\d{3,})(?P<direction>D?)\.nc')
CORE = ('PRES', 'TEMP', 'PSAL')
MEASUREMENT_QC = set('01234589 ')

# Variable: (NetCDF datatype, ordered dimensions). Optional extensions are
# inventoried, but are not treated as mandatory core-profile variables.
SCHEMA = {}
for name, width in {
    'DATA_TYPE': 'STRING16', 'FORMAT_VERSION': 'STRING4', 'HANDBOOK_VERSION': 'STRING4',
    'REFERENCE_DATE_TIME': 'DATE_TIME', 'DATE_CREATION': 'DATE_TIME', 'DATE_UPDATE': 'DATE_TIME',
}.items():
    SCHEMA[name] = ('S1', (width,))
for name, width in {
    'PLATFORM_NUMBER': 'STRING8', 'PROJECT_NAME': 'STRING64', 'PI_NAME': 'STRING64',
    'DATA_CENTRE': 'STRING2', 'DC_REFERENCE': 'STRING32', 'DATA_STATE_INDICATOR': 'STRING4',
    'PLATFORM_TYPE': 'STRING32', 'FLOAT_SERIAL_NO': 'STRING32', 'FIRMWARE_VERSION': 'STRING32',
    'WMO_INST_TYPE': 'STRING4', 'POSITIONING_SYSTEM': 'STRING8',
    'VERTICAL_SAMPLING_SCHEME': 'STRING256',
}.items():
    SCHEMA[name] = ('S1', ('N_PROF', width))
for name in ('DIRECTION', 'DATA_MODE', 'JULD_QC', 'POSITION_QC'):
    SCHEMA[name] = ('S1', ('N_PROF',))
for name in ('CYCLE_NUMBER', 'CONFIG_MISSION_NUMBER'):
    SCHEMA[name] = ('i4', ('N_PROF',))
for name in ('JULD', 'JULD_LOCATION', 'LATITUDE', 'LONGITUDE'):
    SCHEMA[name] = ('f8', ('N_PROF',))
SCHEMA['STATION_PARAMETERS'] = ('S1', ('N_PROF', 'N_PARAM', 'STRING16'))
SCHEMA['PARAMETER'] = ('S1', ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING16'))
for name in ('EQUATION', 'COEFFICIENT', 'COMMENT', 'DATE'):
    SCHEMA['SCIENTIFIC_CALIB_' + name] = (
        'S1', ('N_PROF', 'N_CALIB', 'N_PARAM', 'DATE_TIME' if name == 'DATE' else 'STRING256'))
for name, width in {
    'INSTITUTION': 'STRING4', 'STEP': 'STRING4', 'SOFTWARE': 'STRING4',
    'SOFTWARE_RELEASE': 'STRING4', 'REFERENCE': 'STRING64', 'DATE': 'DATE_TIME',
    'ACTION': 'STRING4', 'PARAMETER': 'STRING16', 'QCTEST': 'STRING16',
}.items():
    SCHEMA['HISTORY_' + name] = ('S1', ('N_HISTORY', 'N_PROF', width))
for name in ('START_PRES', 'STOP_PRES', 'PREVIOUS_VALUE'):
    SCHEMA['HISTORY_' + name] = ('f4', ('N_HISTORY', 'N_PROF'))
for parameter in CORE:
    SCHEMA['PROFILE_' + parameter + '_QC'] = ('S1', ('N_PROF',))
    for suffix in ('', '_ADJUSTED', '_ADJUSTED_ERROR'):
        SCHEMA[parameter + suffix] = ('f4', ('N_PROF', 'N_LEVELS'))
    for suffix in ('_QC', '_ADJUSTED_QC'):
        SCHEMA[parameter + suffix] = ('S1', ('N_PROF', 'N_LEVELS'))


def _text(value):
    return np.asarray(value).tobytes().decode('ascii', errors='replace').strip(' \x00')


def _histogram(values):
    return dict(sorted(Counter(_text(v) or 'blank' for v in np.asarray(values).ravel()).items()))


class FileCheck:
    def __init__(self, path, now=None):
        self.path = Path(path)
        self.now = now or datetime.now(timezone.utc)
        self.result = {'file': self.path.name, 'path': str(self.path), 'profiles': [], 'issues': []}
        self.arrays = {}
        self.good_schema = set()
        self.update = None

    def issue(self, code, message, variable=None, profile=None, mask=None, severity='error'):
        entry = {'severity': severity, 'code': code, 'message': message}
        if variable is not None:
            entry['variable'] = variable
        if profile is not None:
            entry['profile_index'] = int(profile)
        if mask is not None:
            where = np.argwhere(mask)
            if not len(where):
                return
            entry['count'] = len(where)
            entry['sample_indices'] = where[:5].tolist()
        self.result['issues'].append(entry)

    def numeric(self, name, ip=None):
        array = self.arrays[name] if ip is None else self.arrays[name][ip]
        present = np.isfinite(array)
        fill = getattr(self.ds[name], '_FillValue', None)
        if fill is not None:
            present &= array != fill
        return array, present

    def date(self, value, variable, ip=None, required=True, reference=False):
        value = _text(value)
        if not value and not required:
            return None
        try:
            if not re.fullmatch(r'\d{14}', value):
                raise ValueError
            date = datetime.strptime(value, '%Y%m%d%H%M%S').replace(tzinfo=timezone.utc)
        except ValueError:
            self.issue('INVALID_DATE', f'Expected YYYYMMDDHHMISS, got {value!r}', variable, ip)
            return None
        if not reference:
            if date < datetime(1997, 1, 1, tzinfo=timezone.utc) or date > self.now:
                self.issue('DATE_RANGE', f'Date {value} is before 1997 or in the future', variable, ip)
            if self.update is not None and date > self.update:
                self.issue('DATE_AFTER_UPDATE', f'{value} is later than DATE_UPDATE', variable, ip)
        return date

    def schema(self):
        ds = self.ds
        self.result['netcdf_format'] = ds.file_format
        self.result['dimensions'] = {name: len(dim) for name, dim in ds.dimensions.items()}
        self.result['variables'] = {}
        if ds.groups:
            self.issue('UNSUPPORTED_GROUPS', 'Nested NetCDF groups are outside core profile coverage')
        for name in ('title', 'institution', 'source', 'history', 'references', 'Conventions', 'user_manual_version'):
            if name not in ds.ncattrs() or not str(ds.getncattr(name)).strip():
                self.issue('MISSING_GLOBAL_ATTRIBUTE', f'Missing/empty global attribute {name}')
        for name in ('N_PROF', 'N_LEVELS', 'N_PARAM', 'N_CALIB', 'N_HISTORY'):
            if name not in ds.dimensions or len(ds.dimensions[name]) == 0:
                self.issue('MISSING_DIMENSION', f'{name} must exist and be nonempty')
        for name, size in {'DATE_TIME': 14, **{f'STRING{i}': i for i in (2, 4, 8, 16, 32, 64, 256)}}.items():
            if name not in ds.dimensions or len(ds.dimensions[name]) != size:
                self.issue('DIMENSION_SIZE', f'{name} must have length {size}')
        for name, var in ds.variables.items():
            try:
                values = np.asarray(var[:])
                self.arrays[name] = values
            except (ValueError, OSError, RuntimeError) as exc:
                self.issue('UNREADABLE_VARIABLE', str(exc), name)
                continue
            self.result['variables'][name] = {'dtype': str(var.dtype), 'dimensions': list(var.dimensions),
                                               'shape': list(var.shape)}
            if np.issubdtype(values.dtype, np.number):
                self.issue('NONFINITE_VALUE', 'Use the declared fill value instead of NaN/Inf', name,
                           mask=~np.isfinite(values))
            elif values.dtype.kind == 'S':
                self.issue('NULL_CHARACTER', 'Character variables must use spaces, not NULL padding', name,
                           mask=values == b'')
        expected = dict(SCHEMA)
        if 'PARAMETER_DATA_MODE' in ds.variables:
            expected['PARAMETER_DATA_MODE'] = ('S1', ('N_PROF', 'N_PARAM'))
        for name, (dtype, dims) in expected.items():
            if name not in ds.variables:
                self.issue('MISSING_VARIABLE', 'Required core profile variable is missing', name)
                continue
            var = ds[name]
            if var.dimensions != dims or np.dtype(var.dtype).kind != np.dtype(dtype).kind or np.dtype(var.dtype).itemsize != np.dtype(dtype).itemsize:
                self.issue('VARIABLE_SCHEMA', f'Expected {dtype} {dims}, got {var.dtype} {var.dimensions}', name)
            elif name in self.arrays:
                self.good_schema.add(name)
            if name in CORE or name.endswith(('_ADJUSTED', '_ADJUSTED_ERROR')):
                if not hasattr(var, '_FillValue'):
                    self.issue('MISSING_FILL_VALUE', 'Measurement variable requires _FillValue', name)
                for attribute in ('units', 'long_name'):
                    if not hasattr(var, attribute):
                        self.issue('MISSING_VARIABLE_ATTRIBUTE', f'Missing {attribute}', name, severity='warning')
        if 'FORMAT_VERSION' in self.good_schema:
            version = _text(self.arrays['FORMAT_VERSION'])
            self.result['format_version'] = version
            if version != '3.1':
                self.issue('UNSUPPORTED_VERSION', f'Coverage is for profile format 3.1; file declares {version!r}', severity='warning')
        if 'DATA_TYPE' in self.good_schema and _text(self.arrays['DATA_TYPE']) != 'Argo profile':
            self.issue('DATA_TYPE', 'Expected DATA_TYPE = Argo profile', 'DATA_TYPE')

    def metadata(self):
        a = self.arrays
        if 'DATE_UPDATE' in self.good_schema:
            self.update = self.date(a['DATE_UPDATE'], 'DATE_UPDATE')
        if 'DATE_CREATION' in self.good_schema:
            self.date(a['DATE_CREATION'], 'DATE_CREATION')
        reference = None
        if 'REFERENCE_DATE_TIME' in self.good_schema:
            reference = self.date(a['REFERENCE_DATE_TIME'], 'REFERENCE_DATE_TIME', reference=True)
            if _text(a['REFERENCE_DATE_TIME']) != '19500101000000':
                self.issue('REFERENCE_DATE', 'Expected the Argo epoch 19500101000000', 'REFERENCE_DATE_TIME')
        match = FILE_NAME.fullmatch(self.path.name)
        if not match:
            self.issue('FILENAME', 'Expected D<WMO>_<cycle>[D].nc for a core single-cycle file')
        if 'DATA_MODE' in self.good_schema:
            modes = [_text(v) for v in a['DATA_MODE']]
            self.result['data_modes'] = modes
            if 'D' not in modes:
                self.issue('NO_DELAYED_PROFILE', 'A D-file must contain at least one D-mode profile', 'DATA_MODE')
            for ip, mode in enumerate(modes):
                if mode not in ('R', 'A', 'D'):
                    self.issue('DATA_MODE', f'Invalid data mode {mode!r}', 'DATA_MODE', ip)
                elif mode != 'D':
                    self.issue('NON_D_PROFILE', f'This profile is still {mode}-mode; full DM checks skipped', 'DATA_MODE', ip, severity='warning')
        for ip in range(len(self.ds.dimensions.get('N_PROF', ()))):
            is_delayed = 'DATA_MODE' in self.good_schema and _text(a['DATA_MODE'][ip]) == 'D'
            if match and 'PLATFORM_NUMBER' in self.good_schema and _text(a['PLATFORM_NUMBER'][ip]) != match['float']:
                self.issue('PLATFORM_ID', 'Platform number does not match filename', 'PLATFORM_NUMBER', ip)
            if match and 'CYCLE_NUMBER' in self.good_schema:
                value, present = self.numeric('CYCLE_NUMBER', ip)
                if not present or int(value) != int(match['cycle']):
                    self.issue('CYCLE_NUMBER', 'Cycle number does not match filename', 'CYCLE_NUMBER', ip)
            if 'DIRECTION' in self.good_schema:
                direction = _text(a['DIRECTION'][ip])
                if direction not in ('A', 'D') or (match and direction != ('D' if match['direction'] else 'A')):
                    self.issue('DIRECTION', 'Direction is invalid or inconsistent with filename', 'DIRECTION', ip)
            for name in ('JULD_QC', 'POSITION_QC'):
                allowed = set('1234589') if is_delayed else set('01234589 ')
                if name in self.good_schema and _text(a[name][ip]) not in allowed:
                    self.issue('POSITION_TIME_QC', 'Invalid position/time QC (D-mode must be assessed)', name, ip)
            for name, lower, upper in (('LATITUDE', -90, 90), ('LONGITUDE', -180, 180)):
                if name in self.good_schema:
                    value, present = self.numeric(name, ip)
                    if present and not lower <= value <= upper:
                        self.issue('COORDINATE_RANGE', f'Expected {lower} to {upper}', name, ip)
                    if not present and 'POSITION_QC' in self.good_schema and _text(a['POSITION_QC'][ip]) != '9':
                        self.issue('MISSING_POSITION', 'Missing coordinate requires POSITION_QC 9', name, ip)
            for name in ('JULD', 'JULD_LOCATION'):
                if name not in self.good_schema:
                    continue
                value, present = self.numeric(name, ip)
                if not present:
                    if name == 'JULD' and 'JULD_QC' in self.good_schema and _text(a['JULD_QC'][ip]) != '9':
                        self.issue('MISSING_TIME', 'Missing JULD requires JULD_QC 9', name, ip)
                    continue
                if reference is not None:
                    try:
                        date = reference + timedelta(days=float(value))
                        if date < datetime(1997, 1, 1, tzinfo=timezone.utc) or date > self.now:
                            self.issue('DATE_RANGE', 'Observation time is before 1997 or in the future', name, ip)
                        if self.update and date > self.update + timedelta(seconds=1):
                            self.issue('DATE_AFTER_UPDATE', 'Observation time is later than DATE_UPDATE', name, ip)
                    except OverflowError:
                        self.issue('INVALID_DATE', 'Observation time is outside datetime range', name, ip)
            self.profile(ip)

    def profile(self, ip):
        a = self.arrays
        mode = _text(a['DATA_MODE'][ip]) if 'DATA_MODE' in self.good_schema else '?'
        info = {'profile_index': ip, 'data_mode': mode, 'parameters': {}}
        parameter_modes = {}
        if 'VERTICAL_SAMPLING_SCHEME' in self.good_schema:
            info['sampling_scheme'] = _text(a['VERTICAL_SAMPLING_SCHEME'][ip])
        if 'STATION_PARAMETERS' in self.good_schema:
            names = [_text(v) for v in a['STATION_PARAMETERS'][ip] if _text(v)]
            info['station_parameters'] = names
            if len(names) != len(set(names)):
                self.issue('DUPLICATE_PARAMETER', 'Duplicate STATION_PARAMETERS entries', 'STATION_PARAMETERS', ip)
            for p in CORE:
                if p not in names:
                    self.issue('STATION_PARAMETER', f'{p} missing from STATION_PARAMETERS', 'STATION_PARAMETERS', ip)
            extras = set(names) - set(CORE)
            if extras:
                self.issue('EXTRA_PARAMETERS', f'Only PRES/TEMP/PSAL content checks implemented; extra parameters: {sorted(extras)}', severity='warning', profile=ip)
        if mode == 'D' and 'DATA_STATE_INDICATOR' in self.good_schema:
            state = _text(a['DATA_STATE_INDICATOR'][ip])
            if state != '2C':
                self.issue('DATA_STATE', f'Expected 2C for DM data, got {state!r}', 'DATA_STATE_INDICATOR', ip)
        if 'PARAMETER_DATA_MODE' in self.good_schema:
            modes = [_text(v) for v in a['PARAMETER_DATA_MODE'][ip]]
            info['parameter_data_modes'] = modes
            if 'STATION_PARAMETERS' in self.good_schema:
                active = [_text(v) != '' for v in a['STATION_PARAMETERS'][ip]]
                bad = [v for v, used in zip(modes, active) if used and v not in ('R', 'A', 'D')]
                if bad:
                    self.issue('PARAMETER_MODE', f'Invalid parameter data modes: {bad}', 'PARAMETER_DATA_MODE', ip)
                used_modes = [v for v, used in zip(modes, active) if used]
                parameter_modes = {_text(p): m for p, m in zip(a['STATION_PARAMETERS'][ip], modes) if _text(p)}
                expected = next((m for m in ('D', 'A', 'R') if m in used_modes), None)
                if expected is not None and mode != expected:
                    self.issue('MODE_CONSISTENCY', f'DATA_MODE {mode} disagrees with parameter modes', 'PARAMETER_DATA_MODE', ip)
        for p in CORE:
            required = [p, p + '_QC', p + '_ADJUSTED', p + '_ADJUSTED_QC', p + '_ADJUSTED_ERROR', 'PROFILE_' + p + '_QC']
            if not all(name in self.good_schema for name in required):
                continue
            raw, raw_present = self.numeric(p, ip)
            adjusted, adj_present = self.numeric(p + '_ADJUSTED', ip)
            errors, err_present = self.numeric(p + '_ADJUSTED_ERROR', ip)
            q = a[p + '_ADJUSTED_QC'][ip]
            raw_q = a[p + '_QC'][ip]
            info['parameters'][p] = {
                'raw_count': int(raw_present.sum()), 'adjusted_count': int(adj_present.sum()),
                'error_count': int(err_present.sum()), 'adjusted_qc_counts': _histogram(q),
                'raw_qc_counts': _histogram(raw_q), 'profile_qc': _text(a['PROFILE_' + p + '_QC'][ip]),
            }
            for name, flags in ((p + '_QC', raw_q), (p + '_ADJUSTED_QC', q)):
                self.issue('INVALID_QC', 'Invalid Argo measurement QC code', name, ip,
                           mask=~np.isin(flags, list(c.encode() for c in MEASUREMENT_QC)))
            if parameter_modes.get(p, mode) != 'D':
                if mode == 'D':
                    self.issue('NON_D_PARAMETER', f'{p} is not D-mode; DM content checks skipped',
                               'PARAMETER_DATA_MODE', ip, severity='warning')
                continue
            needs_data = np.isin(q, [b'1', b'2', b'3', b'5', b'8'])
            self.issue('UNASSESSED_QC', 'D-mode adjusted QC cannot be 0', p + '_ADJUSTED_QC', ip, mask=q == b'0')
            self.issue('MISSING_QC', 'Present measurements require adjusted QC', p + '_ADJUSTED_QC', ip,
                       mask=(raw_present | adj_present) & np.isin(q, [b' ', b'']))
            self.issue('MISSING_ADJUSTED', 'QC requires an adjusted value', p + '_ADJUSTED', ip, mask=needs_data & ~adj_present)
            self.issue('MISSING_ERROR', 'QC requires an adjusted uncertainty estimate', p + '_ADJUSTED_ERROR', ip, mask=needs_data & ~err_present)
            rejected = np.isin(q, [b'4', b'9'])
            self.issue('BAD_VALUE_NOT_FILLED', 'Adjusted values with QC 4/9 must be fill values', p + '_ADJUSTED', ip, mask=rejected & adj_present)
            self.issue('BAD_ERROR_NOT_FILLED', 'Errors with QC 4/9 must be fill values', p + '_ADJUSTED_ERROR', ip, mask=rejected & err_present)
            self.issue('NEGATIVE_ERROR', 'Uncertainty cannot be negative', p + '_ADJUSTED_ERROR', ip, mask=err_present & (errors < 0))
            self.issue('ORPHAN_ERROR', 'Uncertainty has no corresponding adjusted value', p + '_ADJUSTED_ERROR', ip, mask=err_present & ~adj_present)
            # Independent reference-table calculation; do not reuse writer code.
            assessed = np.isin(q, [b'1', b'2', b'3', b'4', b'5', b'8'])
            good = np.isin(q, [b'1', b'2', b'5', b'8'])
            fraction = good.sum() / assessed.sum() if assessed.any() else None
            expected = ' ' if fraction is None else ('A' if fraction == 1 else 'B' if fraction >= .75 else 'C' if fraction >= .5 else 'D' if fraction >= .25 else 'E' if fraction > 0 else 'F')
            if _text(a['PROFILE_' + p + '_QC'][ip]) != expected.strip():
                self.issue('PROFILE_QC', f'Profile QC should be {expected!r}', 'PROFILE_' + p + '_QC', ip)
        if mode == 'D':
            for p in ('TEMP', 'PSAL'):
                if all(name in self.good_schema for name in ('PRES_ADJUSTED_QC', p + '_ADJUSTED_QC')):
                    bad_pressure = a['PRES_ADJUSTED_QC'][ip] == b'4'
                    usable = np.isin(a[p + '_ADJUSTED_QC'][ip], [b'1', b'2', b'3', b'5', b'8'])
                    self.issue('PRESSURE_QC_DEPENDENCY', 'Bad pressure cannot accompany usable T/S', p + '_ADJUSTED_QC', ip, mask=bad_pressure & usable)
            self.calibration(ip)
        self.result['profiles'].append(info)

    def calibration(self, ip):
        a = self.arrays
        if 'PARAMETER' in self.good_schema:
            recorded = set()
            for record, row in enumerate(a['PARAMETER'][ip]):
                for slot, value in enumerate(row):
                    p = _text(value)
                    if not p:
                        continue
                    recorded.add(p)
                    for field in ('EQUATION', 'COEFFICIENT', 'COMMENT', 'DATE'):
                        name = 'SCIENTIFIC_CALIB_' + field
                        if name not in self.good_schema:
                            continue
                        cell = a[name][ip, record, slot]
                        if not _text(cell):
                            self.issue('MISSING_CALIBRATION', f'{p}, calibration {record}: {field} is empty', name, ip)
                        elif field == 'DATE':
                            self.date(cell, name, ip)
            if 'STATION_PARAMETERS' in self.good_schema:
                for cell in a['STATION_PARAMETERS'][ip]:
                    p = _text(cell)
                    if p and p not in recorded:
                        self.issue('MISSING_CALIBRATION_PARAMETER', f'No calibration record for {p}', 'PARAMETER', ip)
        if 'HISTORY_DATE' in self.good_schema:
            populated = False
            dm_record = False
            for record, cell in enumerate(a['HISTORY_DATE'][:, ip]):
                if _text(cell):
                    populated = True
                    self.date(cell, 'HISTORY_DATE', ip)
                    if 'HISTORY_STEP' in self.good_schema and _text(a['HISTORY_STEP'][record, ip]) == 'ARSQ':
                        dm_record = True
            if not populated:
                self.issue('MISSING_HISTORY', 'No populated history record for this profile', 'HISTORY_DATE', ip)
            elif not dm_record:
                self.issue('MISSING_DM_HISTORY', 'No ARSQ delayed-mode history entry', 'HISTORY_STEP', ip, severity='warning')

    def run(self):
        try:
            with Dataset(self.path, 'r') as ds:
                self.ds = ds
                ds.set_auto_maskandscale(False)
                ds.set_auto_chartostring(False)
                self.schema()
                self.metadata()
        except (OSError, RuntimeError, ValueError, IndexError, KeyError, TypeError, OverflowError) as exc:
            self.issue('FILE_CHECK_FAILED', f'Could not finish inspecting this file: {type(exc).__name__}: {exc}')
        errors = sum(i['severity'] == 'error' for i in self.result['issues'])
        warnings = sum(i['severity'] == 'warning' for i in self.result['issues'])
        self.result.update(error_count=errors, warning_count=warnings,
                           status='FAIL' if errors else 'WARN' if warnings else 'PASS')
        return self.result


def verify_file(path):
    return FileCheck(path).run()


def verify_directory(path):
    """Inspect a float folder or D/ folder; never open any file for writing."""
    path = Path(path)
    if not path.is_dir():
        raise ValueError(f'Not a directory: {path}')
    directory = path / 'D' if (path / 'D').is_dir() else path
    files = sorted(directory.glob('D*.nc'))
    if not files:
        raise ValueError(f'No D*.nc files found in {directory}')
    results = [verify_file(file) for file in files]
    directory_issues = []
    source_dir = directory.parent / 'R'
    if directory.name == 'D' and source_dir.is_dir():
        missing = [f'D{p.name[1:]}' for p in sorted(source_dir.glob('R*.nc'))
                   if FILE_NAME.fullmatch(f'D{p.name[1:]}') and not (directory / f'D{p.name[1:]}').is_file()]
        if missing:
            directory_issues.append({'severity': 'warning', 'code': 'MISSING_D_FILES',
                                     'message': 'Some local R-files have no D-file; this may be an intentional subset',
                                     'files': missing})
    codes = Counter(i['code'] for result in results for i in result['issues'])
    codes.update(i['code'] for i in directory_issues)
    return {
        'schema_version': 1, 'checker': 'verify_dfiles', 'scope': SCOPE, 'sources': SOURCES,
        'created_utc': datetime.now(timezone.utc).isoformat(), 'directory': str(directory),
        'summary': {'files': len(results), 'profiles': sum(len(r['profiles']) for r in results),
                    'passed': sum(r['status'] == 'PASS' for r in results),
                    'failed': sum(r['status'] == 'FAIL' for r in results),
                    'warning_only': sum(r['status'] == 'WARN' for r in results),
                    'errors': sum(r['error_count'] for r in results),
                    'warnings': sum(r['warning_count'] for r in results) + len(directory_issues),
                    'issue_counts': dict(sorted(codes.items()))},
        'directory_issues': directory_issues, 'files': results,
    }
