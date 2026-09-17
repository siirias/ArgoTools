"""Write all resolved profiles, applying rejections and retaining accepted data.

Raw measurements/QC remain untouched. Rejected adjusted data and errors are
filled, as required for Argo adjusted QC 4; no numerical correction is computed.
"""
from collections import defaultdict
from datetime import datetime, timezone
import yaml
import os
from pathlib import Path
import tempfile

import numpy as np
from netCDF4 import Dataset

from .download import file_identity
from .profiles import core_parameters
from .instructions import CORE_PARAMETERS, sha256, uncertainty_value


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
        raise ValueError('Only whole-profile PRES/TEMP/PSAL processing is implemented')


def validate_source(ds, decision):
    target = decision.target
    validate_target(ds, target)
    for suggestion in decision.suggestions:
        if type(suggestion.priority) is not int:
            raise ValueError("priority must be an integer")
        st = suggestion.target
        if (st.source != target.source or st.profile_index != target.profile_index or
                st.selection != 'whole_profile' or st.sample_indices is not None or st.pressure_range is not None):
            raise ValueError('Decision contains inconsistent targets')
        if suggestion.action == 'set_uncertainty':
            uncertainty_value(suggestion.value)
            if not st.parameters or any(p not in CORE_PARAMETERS for p in st.parameters) or suggestion.flag is not None:
                raise ValueError('Unsupported uncertainty instruction')
        elif (suggestion.action not in ('flag', 'accept', 'no_finding') or st.parameters != CORE_PARAMETERS or
              (suggestion.action == 'flag' and suggestion.flag != '4')):
            raise ValueError('Decision contains unsupported instructions')
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
    active = core_parameters(ds, target.profile_index)
    if not active:
        raise ValueError(f'{target.source} profile {target.profile_index}: no supported core parameters')
    for p in CORE_PARAMETERS:
        if p not in active:
            for suffix in ('', '_ADJUSTED', '_ADJUSTED_ERROR'):
                if p + suffix in ds.variables and valid(ds[p + suffix], target.profile_index).any():
                    raise ValueError(f'{target.source} profile {target.profile_index}: {p} data exist but are not declared')
    for p in active:
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
    for p in active:
        value = decision.uncertainty(p)
        if value is not None:
            var = ds[p + '_ADJUSTED_ERROR']
            if value >= np.finfo(var.dtype).max or value < np.finfo(var.dtype).tiny or np.asarray(value, dtype=var.dtype) == var._FillValue:
                raise ValueError(f'{p}: uncertainty is not representable as a non-fill value')
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


def retained_arrays(ds, parameter, iprof, uncertainty=None):
    """Retain existing adjustments/QC; use raw data when adjustments are absent.

    Existing flags 3/4/9 are not upgraded. Fill values for QC 4/9 are enforced;
    no corrections or uncertainty estimates are calculated.
    """
    adjusted = np.ma.masked_invalid(np.ma.asarray(ds[parameter + '_ADJUSTED'][iprof]).copy())
    errors = np.ma.masked_invalid(np.ma.asarray(ds[parameter + '_ADJUSTED_ERROR'][iprof]).copy())
    raw = np.ma.masked_invalid(np.ma.asarray(ds[parameter][iprof]))
    flags = np.ma.filled(ds[parameter + '_ADJUSTED_QC'][iprof], b' ').copy()
    raw_flags = np.ma.filled(ds[parameter + '_QC'][iprof], b'0')
    unspecified = np.isin(flags, [b' ', b'0'])
    copied = np.ma.getmaskarray(adjusted) & ~np.ma.getmaskarray(raw) & unspecified
    adjusted[copied] = raw[copied]
    flags[copied | (flags == b' ')] = raw_flags[copied | (flags == b' ')]
    flags[np.ma.getmaskarray(adjusted) & (flags != b'4')] = b'9'
    flags[flags == b' '] = b'0'
    missing_or_bad = np.isin(flags, [b'4', b'9'])
    adjusted[missing_or_bad] = np.ma.masked
    errors[missing_or_bad | np.ma.getmaskarray(adjusted) | copied] = np.ma.masked
    if uncertainty is not None:
        eligible = ~np.ma.getmaskarray(adjusted) & np.isin(flags, [b'1', b'2', b'3', b'5', b'8'])
        errors[eligible] = uncertainty
    return adjusted, flags, errors, int(copied.sum())


def profile_qc(flags):
    """Argo profile summary: fraction of nonmissing levels with QC 1/2/5/8."""
    count = int(np.count_nonzero(flags != b'9'))
    if count == 0:
        return b' '
    good = int(np.count_nonzero(np.isin(flags, [b'1', b'2', b'5', b'8'])))
    if good == count:
        return b'A'
    if good == 0:
        return b'F'
    fraction = good / count
    return b'B' if fraction >= .75 else b'C' if fraction >= .5 else b'D' if fraction >= .25 else b'E'


def previous_calibration(ds, parameter, iprof, before):
    """Keep the last recorded adjustment formula when retaining adjusted data."""
    for record in reversed(range(before)):
        for slot, name in enumerate(ds['PARAMETER'][iprof, record]):
            if text(name) == parameter:
                equation = text(ds['SCIENTIFIC_CALIB_EQUATION'][iprof, record, slot])
                coefficient = text(ds['SCIENTIFIC_CALIB_COEFFICIENT'][iprof, record, slot])
                if equation or coefficient:
                    return equation or 'none', coefficient or 'none'
    return 'none', 'none'


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
        reasons = '; '.join(f'{s.checker}: {s.reason}' for s in decision.winning_suggestions if s.action == 'flag')
        for p in core_parameters(ds, ip):
            previous = np.ma.filled(ds[p + '_ADJUSTED_QC'][ip], b' ')
            copied = 0
            if decision.rejected:
                present = observed_samples(ds, p, ip)
                flags = np.where(present, b'4', b'9')
                for suffix in ('_ADJUSTED', '_ADJUSTED_ERROR'):
                    var = ds[p + suffix]
                    var[ip, :] = np.ma.masked_all(var.shape[1], dtype=var.dtype)
                comment = f'Bad profile; no correction applied. {reasons}'
            else:
                adjusted, flags, errors, copied = retained_arrays(ds, p, ip, decision.uncertainty(p))
                ds[p + '_ADJUSTED'][ip, :] = adjusted
                ds[p + '_ADJUSTED_ERROR'][ip, :] = errors
                comment = ('Profile retained by resolved instructions. Existing adjustments and QC retained; '
                           'raw data used where adjustments absent.')
                accepted = [s for s in decision.winning_suggestions if s.action == 'accept']
                if accepted:
                    comment = f'Accepted at priority {accepted[0].priority}. ' + comment
                estimate = decision.uncertainty(p)
                if estimate is not None:
                    scope = 'profile' if any(s.action == 'set_uncertainty' and p in s.target.parameters
                                             and s.scope == 'profile' for s in decision.suggestions) else 'float'
                    provenance = '; '.join(f'{s.checker}: {s.reason}' for s in decision.suggestions
                                           if s.action == 'set_uncertainty' and p in s.target.parameters and s.scope == scope)
                    comment = f'Assigned {p}_ADJUSTED_ERROR={estimate:g}. {provenance}. ' + comment
                else:
                    comment += ' No uncertainty estimated.'
            ds[p + '_ADJUSTED_QC'][ip, :] = flags
            ds['PROFILE_' + p + '_QC'][ip] = profile_qc(flags)
            j = names.index(p)
            equation, coefficient = ('none', 'none') if decision.rejected else previous_calibration(ds, p, ip, calib_index)
            if copied:
                equation = f'{p}_ADJUSTED = {p} for samples copied from raw; existing adjustments retained elsewhere'
                coefficient = 'none'
            if 'PARAMETER_DATA_MODE' in ds.variables:
                ds['PARAMETER_DATA_MODE'][ip, j] = b'D'
            for name, value in {
                'PARAMETER': p,
                'SCIENTIFIC_CALIB_DATE': now,
                'SCIENTIFIC_CALIB_EQUATION': equation,
                'SCIENTIFIC_CALIB_COEFFICIENT': coefficient,
                'SCIENTIFIC_CALIB_COMMENT': comment,
            }.items():
                put_text(ds[name], (ip, calib_index, j), value)
            fields = {
                'DATE': now, 'ACTION': 'CF' if (previous != flags).any() else '',
                'STEP': 'ARSQ', 'PARAMETER': p,
                'INSTITUTION': meta['institution'], 'SOFTWARE': meta['software'],
                'SOFTWARE_RELEASE': meta['software_release'], 'REFERENCE': report_name,
            }
            for field, value in fields.items():
                put_text(ds['HISTORY_' + field], (history, ip), value)
            if 'HISTORY_COMMENT' in ds.variables:
                require(ds, 'HISTORY_COMMENT', ('N_HISTORY', 'N_PROF'), char=True)
                put_text(ds['HISTORY_COMMENT'], (history, ip), comment)
            history += 1
            changes.append({'profile_index': ip, 'parameter': p,
                            'outcome': 'reject' if decision.rejected else 'retain',
                            'bad_samples': int((flags == b'4').sum()),
                            'missing_samples': int((flags == b'9').sum()),
                            'raw_samples_copied': copied,
                            'samples_without_uncertainty': int((valid(ds[p + '_ADJUSTED'], ip) &
                                                                ~valid(ds[p + '_ADJUSTED_ERROR'], ip)).sum()),
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
            for p in core_parameters(src, ip):
                if decision.rejected:
                    expected = np.where(observed_samples(src, p, ip), b'4', b'9')
                    if valid(dst[p + '_ADJUSTED'], ip).any() or valid(dst[p + '_ADJUSTED_ERROR'], ip).any():
                        raise ValueError('Rejected adjusted values/errors must be missing')
                else:
                    adjusted, expected, errors, _ = retained_arrays(src, p, ip, decision.uncertainty(p))
                    for suffix, values in (('_ADJUSTED', adjusted), ('_ADJUSTED_ERROR', errors)):
                        actual = dst[p + suffix][ip]
                        np.testing.assert_array_equal(np.ma.getmaskarray(actual), np.ma.getmaskarray(values))
                        np.testing.assert_array_equal(np.ma.filled(actual, np.nan), np.ma.filled(values, np.nan))
                np.testing.assert_array_equal(dst[p + '_ADJUSTED_QC'][ip], expected)
                np.testing.assert_array_equal(dst['PROFILE_' + p + '_QC'][ip], profile_qc(expected))
        active_by_profile = {ip: core_parameters(src, ip) for ip in profiles}
        # Compare stored data, including absent parameters and existing metadata records.
        src.set_auto_maskandscale(False)
        dst.set_auto_maskandscale(False)
        for name, var in src.variables.items():
            if name == 'DATE_UPDATE':
                continue
            old = var[...]
            new = dst[name][tuple(slice(0, n) for n in var.shape)]
            if name in profile_edits:
                parameter = next((p for p in CORE_PARAMETERS if name in
                                  (p + '_ADJUSTED', p + '_ADJUSTED_ERROR', p + '_ADJUSTED_QC', 'PROFILE_' + p + '_QC')), None)
                untouched = [i for i in range(len(src.dimensions['N_PROF']))
                             if i not in profiles or (parameter is not None and parameter not in active_by_profile[i])]
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
    Audit reports in the sibling reports/ directory include hashes and reasons.
    """
    metadata = {**DEFAULT_META, **(meta or {})}
    plan = prepare(decisions, r_dir, d_dir, overwrite)
    reports_dir = Path(d_dir).parent / 'reports'
    for _, output, _, _, _ in plan:
        report_path = reports_dir / output.with_suffix('.report.yaml').name
        if report_path.exists() and not overwrite:
            raise FileExistsError(f'{report_path} exists; use --overwrite')
    if dry_run:
        return [{'source': str(s), 'output': str(o), 'status': 'would_write',
                 'report': str(reports_dir / o.with_suffix('.report.yaml').name),
                 'decisions': [d.as_dict() for d in entries]} for s, o, entries, _, _ in plan]
    if not plan:
        return []
    Path(d_dir).mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    with tempfile.TemporaryDirectory(dir=d_dir, prefix='.dmqc-') as staging:
        staged = []
        for source, output, entries, history, calib in plan:
            report_path = reports_dir / output.with_suffix('.report.yaml').name
            temporary = Path(staging) / output.name
            with Dataset(source) as src:
                history_count = sum(len(core_parameters(src, d.target.profile_index)) for d in entries)
            clone_expanded(source, temporary, history_count)
            with Dataset(temporary, 'r+') as ds:
                ds.set_auto_chartostring(False)
                changes = _apply(ds, entries, metadata, history, calib, f'../reports/{report_path.name}')
            validate_output(source, temporary, entries)
            if sha256(source) != entries[0].source_sha256:
                raise ValueError(f'{source}: source changed during writing')
            report = {'source': str(source), 'output': str(output), 'report': str(report_path), 'status': 'written',
                      'source_sha256': entries[0].source_sha256, 'output_sha256': sha256(temporary),
                      'created_utc': datetime.now(timezone.utc).isoformat(),
                      'metadata': metadata, 'changes': changes,
                      'decisions': [d.as_dict() for d in entries]}
            audit = Path(staging) / report_path.name
            audit.write_text(yaml.safe_dump(report, sort_keys=False, allow_unicode=True), encoding='utf-8')
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
