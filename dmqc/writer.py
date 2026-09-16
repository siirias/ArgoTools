"""Apply resolved whole-profile rejections and atomically publish D-files.

Raw measurements/QC remain untouched. Rejected adjusted data and errors are
filled, as required for Argo adjusted QC 4; no numerical correction is computed.
"""
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from netCDF4 import Dataset

from .download import file_identity
from .instructions import CORE_PARAMETERS, sha256


DEFAULT_META = {'institution': 'IF', 'operator': 'unknown', 'software': 'ADMW',
                'software_release': '0.1'}


def text(value):
    return np.ma.filled(value, b' ').tobytes().decode('ascii', errors='replace').strip(' \x00')


def put_text(variable, index, value):
    length = variable.shape[-1]
    value = str(value).encode('ascii', errors='replace')[:length].ljust(length, b' ')
    variable[index] = np.frombuffer(value, dtype='S1')


def valid(variable, profile):
    values = np.ma.asarray(variable[profile, :])
    return ~np.ma.getmaskarray(values) & np.isfinite(np.ma.filled(values, np.nan))


def require(ds, name, dimensions, char=False):
    if name not in ds.variables:
        raise ValueError(f'Missing required variable {name}')
    v = ds[name]
    if char:
        correct = v.dimensions[:-1] == dimensions and v.dtype == np.dtype('S1')
    else:
        correct = v.dimensions == dimensions
    if not correct:
        raise ValueError(f'{name}: unsupported dimensions/type {v.dimensions}, {v.dtype}')
    return v


def validate_target(ds, target):
    float_id, cycle = file_identity(target.source)
    ip = target.profile_index
    if type(ip) is not int or ip < 0 or ip >= len(ds.dimensions['N_PROF']):
        raise ValueError(f'{target.source}: invalid profile index {ip}')
    platform = require(ds, 'PLATFORM_NUMBER', ('N_PROF',), char=True)
    cycles = require(ds, 'CYCLE_NUMBER', ('N_PROF',))
    if text(platform[ip]) != float_id or np.ma.is_masked(cycles[ip]) or int(cycles[ip]) != cycle:
        raise ValueError(f'{target.source}: filename does not match profile identity')
    if target.selection != 'whole_profile' or target.parameters != CORE_PARAMETERS or target.sample_indices is not None or target.pressure_range is not None:
        raise ValueError('Only whole-profile PRES/TEMP/PSAL rejection is implemented')


def validate_source(ds, decision):
    target = decision.target
    validate_target(ds, target)
    if not decision.suggestions or not any(s.action == 'flag' and s.flag == '4' for s in decision.suggestions):
        raise ValueError('A rejection requires an explicit bad-data instruction')
    if any(s.target != target or s.action not in ('flag', 'no_finding') or
           (s.action == 'flag' and s.flag != '4') for s in decision.suggestions):
        raise ValueError('Decision contains unsupported or inconsistent instructions')
    require(ds, 'DATA_MODE', (), char=True)  # One character per profile.
    if ds['DATA_MODE'].dimensions != ('N_PROF',):
        raise ValueError('DATA_MODE must use N_PROF')
    require(ds, 'DATA_STATE_INDICATOR', ('N_PROF',), char=True)
    require(ds, 'DATE_UPDATE', (), char=True)
    require(ds, 'STATION_PARAMETERS', ('N_PROF', 'N_PARAM'), char=True)
    if 'PARAMETER_DATA_MODE' in ds.variables:
        require(ds, 'PARAMETER_DATA_MODE', ('N_PROF',), char=True)
        if ds['PARAMETER_DATA_MODE'].dimensions != ('N_PROF', 'N_PARAM'):
            raise ValueError('PARAMETER_DATA_MODE must use N_PROF, N_PARAM')
    names = [text(row) for row in ds['STATION_PARAMETERS'][target.profile_index]]
    if any(names.count(p) != 1 for p in CORE_PARAMETERS):
        raise ValueError('Each core parameter must appear once in STATION_PARAMETERS')
    for p in CORE_PARAMETERS:
        for suffix in ('', '_ADJUSTED', '_ADJUSTED_ERROR'):
            v = require(ds, p + suffix, ('N_PROF', 'N_LEVELS'))
            if not np.issubdtype(v.dtype, np.floating) or not hasattr(v, '_FillValue'):
                raise ValueError(f'{p + suffix}: expected floating data with _FillValue')
        require(ds, p + '_QC', ('N_PROF',), char=True)
        require(ds, p + '_ADJUSTED_QC', ('N_PROF',), char=True)
        if ds[p + '_ADJUSTED_QC'].dimensions != ('N_PROF', 'N_LEVELS') or ds[p + '_QC'].dimensions != ('N_PROF', 'N_LEVELS'):
            raise ValueError(f'{p}: QC arrays must use N_PROF, N_LEVELS')
        require(ds, 'PROFILE_' + p + '_QC', (), char=True)
        if ds['PROFILE_' + p + '_QC'].dimensions != ('N_PROF',):
            raise ValueError(f'{p}: profile QC must use N_PROF')
    require(ds, 'PARAMETER', ('N_PROF', 'N_CALIB', 'N_PARAM'), char=True)
    for field in ('DATE', 'EQUATION', 'COEFFICIENT', 'COMMENT'):
        require(ds, 'SCIENTIFIC_CALIB_' + field, ('N_PROF', 'N_CALIB', 'N_PARAM'), char=True)
    for field in ('DATE', 'ACTION', 'STEP', 'PARAMETER', 'INSTITUTION', 'SOFTWARE', 'SOFTWARE_RELEASE', 'REFERENCE'):
        require(ds, 'HISTORY_' + field, ('N_HISTORY', 'N_PROF'), char=True)


def clone_expanded(source, destination, history_count):
    """Preserve stored values/attributes while extending fixed calibration dims."""
    with Dataset(source) as src, Dataset(destination, 'w', format=src.file_format) as dst:
        src.set_auto_maskandscale(False)
        src.set_auto_chartostring(False)
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
        for name, dim in src.dimensions.items():
            size = len(dim)
            if name == 'N_CALIB':
                size += 1
            elif name == 'N_HISTORY':
                size += history_count
            dst.createDimension(name, None if dim.isunlimited() else size)
        for name, var in src.variables.items():
            kwargs = {}
            if '_FillValue' in var.ncattrs():
                kwargs['fill_value'] = var.getncattr('_FillValue')
            copied = dst.createVariable(name, var.datatype, var.dimensions, **kwargs)
            copied.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != '_FillValue'})
            copied.set_auto_maskandscale(False)
            copied.set_auto_chartostring(False)
            copied[tuple(slice(0, n) for n in var.shape)] = var[...]


def observed_samples(ds, parameter, profile):
    # R-mode files may have no adjusted values yet. Raw data identify the samples
    # that were rejected; absent measurements stay missing, not bad.
    present = valid(ds[parameter], profile) | valid(ds[parameter + '_ADJUSTED'], profile)
    raw_qc = np.ma.filled(ds[parameter + '_QC'][profile], b' ')
    adjusted_qc = np.ma.filled(ds[parameter + '_ADJUSTED_QC'][profile], b' ')
    return present & (raw_qc != b'9') & (adjusted_qc != b'9')


def _apply(ds, decisions, meta, history_start, calib_index, report_name):
    now = datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')
    put_text(ds['DATE_UPDATE'], slice(None), now)
    if 'DATE_UPDATE' in ds.ncattrs():
        ds.setncattr('DATE_UPDATE', now)
    changes = []
    history = history_start
    for decision in decisions:
        ip = decision.target.profile_index
        ds['DATA_MODE'][ip] = b'D'
        put_text(ds['DATA_STATE_INDICATOR'], ip, '2C')
        names = [text(row) for row in ds['STATION_PARAMETERS'][ip]]
        reasons = '; '.join(f'{s.checker}: {s.reason}' for s in decision.suggestions if s.action == 'flag')
        for p in CORE_PARAMETERS:
            present = observed_samples(ds, p, ip)
            previous = np.ma.filled(ds[p + '_ADJUSTED_QC'][ip], b' ')
            flags = np.full(present.shape, b'9', dtype='S1')
            flags[present] = b'4'
            ds[p + '_ADJUSTED_QC'][ip, :] = flags
            for suffix in ('_ADJUSTED', '_ADJUSTED_ERROR'):
                var = ds[p + suffix]
                var[ip, :] = np.ma.masked_all(var.shape[1], dtype=var.dtype)
            ds['PROFILE_' + p + '_QC'][ip] = b'F' if present.any() else b' '
            j = names.index(p)
            if 'PARAMETER_DATA_MODE' in ds.variables:
                ds['PARAMETER_DATA_MODE'][ip, j] = b'D'
            for name, value in {
                'PARAMETER': p,
                'SCIENTIFIC_CALIB_DATE': now,
                'SCIENTIFIC_CALIB_EQUATION': 'none',
                'SCIENTIFIC_CALIB_COEFFICIENT': 'none',
                'SCIENTIFIC_CALIB_COMMENT': f'Bad profile; no correction applied. {reasons}',
            }.items():
                put_text(ds[name], (ip, calib_index, j), value)
            fields = {
                'DATE': now, 'ACTION': 'CF', 'STEP': 'ARSQ', 'PARAMETER': p,
                'INSTITUTION': meta['institution'], 'SOFTWARE': meta['software'],
                'SOFTWARE_RELEASE': meta['software_release'], 'REFERENCE': report_name,
            }
            for field, value in fields.items():
                put_text(ds['HISTORY_' + field], (history, ip), value)
            if 'HISTORY_COMMENT' in ds.variables:
                require(ds, 'HISTORY_COMMENT', ('N_HISTORY', 'N_PROF'), char=True)
                put_text(ds['HISTORY_COMMENT'], (history, ip), reasons)
            history += 1
            changes.append({'profile_index': ip, 'parameter': p,
                            'bad_samples': int(present.sum()),
                            'missing_samples': int((~present).sum()),
                            'qc_changes': int((previous != flags).sum())})
    return changes


def validate_output(source, output, decisions):
    """Check intended edits and preserve every original variable outside them."""
    profiles = {d.target.profile_index for d in decisions}
    profile_edits = {'DATA_MODE', 'DATA_STATE_INDICATOR', 'PARAMETER_DATA_MODE'}
    for p in CORE_PARAMETERS:
        profile_edits.update({p + '_ADJUSTED', p + '_ADJUSTED_ERROR', p + '_ADJUSTED_QC', 'PROFILE_' + p + '_QC'})
    with Dataset(source) as src, Dataset(output) as dst:
        src.set_auto_chartostring(False)
        dst.set_auto_chartostring(False)
        for decision in decisions:
            ip = decision.target.profile_index
            if dst['DATA_MODE'][ip] != b'D' or text(dst['DATA_STATE_INDICATOR'][ip]) != '2C':
                raise ValueError('Output mode metadata validation failed')
            for p in CORE_PARAMETERS:
                expected = np.where(observed_samples(src, p, ip), b'4', b'9')
                np.testing.assert_array_equal(dst[p + '_ADJUSTED_QC'][ip], expected)
                if valid(dst[p + '_ADJUSTED'], ip).any() or valid(dst[p + '_ADJUSTED_ERROR'], ip).any():
                    raise ValueError('Rejected adjusted values/errors must be missing')
        # Compare stored data, including fill values and existing metadata records.
        src.set_auto_maskandscale(False)
        dst.set_auto_maskandscale(False)
        for name, var in src.variables.items():
            if name == 'DATE_UPDATE':
                continue
            old = var[...]
            new = dst[name][tuple(slice(0, n) for n in var.shape)]
            if name in profile_edits:
                untouched = [i for i in range(len(src.dimensions['N_PROF'])) if i not in profiles]
                old, new = old[untouched], new[untouched]
            np.testing.assert_array_equal(old, new, err_msg=f'Unexpected change in {name}')


def prepare(decisions, r_dir, d_dir, overwrite=False):
    """Preflight every file before publishing any outputs."""
    grouped = defaultdict(list)
    for decision in decisions:
        grouped[decision.target.source].append(decision)
    plan = []
    for name, entries in sorted(grouped.items()):
        source = Path(r_dir) / name
        output = Path(d_dir) / ('D' + name[1:])
        if output.exists() and not overwrite:
            raise FileExistsError(f'{output} exists; use --overwrite to replace it')
        if source.resolve() == output.resolve():
            raise ValueError('Source and output must differ')
        checksum = sha256(source)
        if any(d.source_sha256 != checksum for d in entries):
            raise ValueError(f'{source}: source changed after combining instructions')
        if len({d.target.profile_index for d in entries}) != len(entries):
            raise ValueError(f'{source}: duplicate resolved profile decisions')
        with Dataset(source) as ds:
            ds.set_auto_chartostring(False)
            for decision in entries:
                validate_source(ds, decision)
            plan.append((source, output, entries, len(ds.dimensions['N_HISTORY']), len(ds.dimensions['N_CALIB'])))
    return plan


def write_d_files(decisions, r_dir, d_dir, meta=None, dry_run=False, overwrite=False):
    """Validate and stage the complete batch, then publish each finished file.

    Publication is atomic per file, not a filesystem-wide batch transaction.
    An audit report alongside each D-file includes source hashes and all reasons.
    """
    metadata = {**DEFAULT_META, **(meta or {})}
    plan = prepare(decisions, r_dir, d_dir, overwrite)
    if dry_run:
        return [{'source': str(s), 'output': str(o), 'status': 'would_write',
                 'decisions': [d.as_dict() for d in entries]} for s, o, entries, _, _ in plan]
    if not plan:
        return []
    Path(d_dir).mkdir(parents=True, exist_ok=True)
    reports = []
    with tempfile.TemporaryDirectory(dir=d_dir, prefix='.dmqc-') as staging:
        staged = []
        for source, output, entries, history, calib in plan:
            report_path = output.with_suffix('.report.json')
            if report_path.exists() and not overwrite:
                raise FileExistsError(f'{report_path} exists; use --overwrite')
            temporary = Path(staging) / output.name
            clone_expanded(source, temporary, len(entries) * len(CORE_PARAMETERS))
            with Dataset(temporary, 'r+') as ds:
                ds.set_auto_chartostring(False)
                changes = _apply(ds, entries, metadata, history, calib, report_path.name)
            validate_output(source, temporary, entries)
            if sha256(source) != entries[0].source_sha256:
                raise ValueError(f'{source}: source changed during writing')
            report = {'source': str(source), 'output': str(output), 'status': 'written',
                      'source_sha256': entries[0].source_sha256, 'output_sha256': sha256(temporary),
                      'created_utc': datetime.now(timezone.utc).isoformat(),
                      'metadata': metadata, 'changes': changes,
                      'decisions': [d.as_dict() for d in entries]}
            audit = Path(staging) / report_path.name
            audit.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
            staged.append((temporary, output, audit, report_path))
            reports.append(report)
        for temporary, output, audit, report_path in staged:
            if overwrite:
                os.replace(temporary, output)
                os.replace(audit, report_path)
            else:
                os.link(audit, report_path)
                try:
                    os.link(temporary, output)
                except BaseException:
                    report_path.unlink()
                    raise
    return reports
