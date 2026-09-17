"""Resolve partner suggestions without performing scientific checks."""
from pathlib import Path
from dataclasses import replace

from .download import file_identity
from .instructions import CORE_PARAMETERS, Decision, Target, FloatTarget, sha256, uncertainty_value


def combine_instructions(instructions, r_dir, float_id=None, cycles=None):
    """Plan every selected R-file/profile; bad wins, others retain existing QC.

    No rejection instruction means retain, not an upgrade to QC 1. Optional
    float/cycle filters restrict both the input instructions and output set.
    """
    from netCDF4 import Dataset
    from .writer import validate_target

    grouped = {}
    hashes = {}
    defaults = []
    for item in instructions:
        if type(item.priority) is not int:
            raise ValueError("priority must be an integer")
        if isinstance(item.target, FloatTarget):
            if item.action != 'set_uncertainty' or item.scope != 'float':
                raise ValueError('Unsupported float-wide operation')
            uncertainty_value(item.value)
            if float_id is not None and item.target.float_id != str(float_id):
                raise ValueError('Float-wide instruction belongs to another float')
            defaults.append(item)
            continue
        target = Target(item.target.source, item.target.profile_index)
        if item.action == 'set_uncertainty':
            uncertainty_value(item.value)
        elif item.target.parameters != CORE_PARAMETERS:
            raise ValueError('Flag decisions must target all core parameters')
        if not item.target.parameters or any(p not in CORE_PARAMETERS for p in item.target.parameters):
            raise ValueError('Unsupported parameters')
        if item.target.selection != 'whole_profile' or item.target.sample_indices is not None or item.target.pressure_range is not None:
            raise ValueError('Only whole-profile core-parameter decisions are implemented')
        if item.action not in ('flag', 'accept', 'no_finding', 'set_uncertainty') or (item.action == 'flag' and item.flag != '4'):
            raise ValueError('Unsupported instruction operation')
        source_float, cycle = file_identity(target.source)
        if float_id is not None and source_float != str(float_id):
            raise ValueError(f'{target.source}: instruction belongs to another float')
        if cycles is not None and cycle not in cycles:
            continue
        if type(target.profile_index) is not int or target.profile_index < 0:
            raise ValueError('profile_index must be a zero-based nonnegative integer')
        if target.source not in hashes:
            hashes[target.source] = sha256(Path(r_dir) / target.source)
        if item.source_sha256 is not None and item.source_sha256 != hashes[target.source]:
            raise ValueError(f'{item.origin}: source changed since the check ran: {target.source}')
        if target not in grouped:
            with Dataset(Path(r_dir) / target.source) as ds:
                ds.set_auto_chartostring(False)
                validate_target(ds, target)
        grouped.setdefault(target, []).append(item)
    # Include every profile, even when the inspector's rejection list is empty.
    decisions = []
    for path in sorted(Path(r_dir).glob('R*.nc')):
        source_float, cycle = file_identity(path.name)
        if float_id is not None and source_float != str(float_id):
            continue
        if cycles is not None and cycle not in cycles:
            continue
        if path.name not in hashes:
            hashes[path.name] = sha256(path)
        with Dataset(path) as ds:
            ds.set_auto_chartostring(False)
            for iprof in range(len(ds.dimensions['N_PROF'])):
                target = Target(path.name, iprof)
                validate_target(ds, target)
                expanded = [replace(item, target=Target(path.name, iprof, item.target.parameters))
                            for item in defaults if item.target.float_id == source_float]
                decision = Decision(target, hashes[path.name], tuple(expanded + grouped.get(target, [])))
                for parameter in CORE_PARAMETERS:
                    decision.uncertainty(parameter)  # Conflicts stop before any output is written.
                decisions.append(decision)
    return decisions
