"""Resolve partner suggestions without performing scientific checks."""
from pathlib import Path

from .download import file_identity
from .instructions import CORE_PARAMETERS, Decision, Target, sha256


def combine_instructions(instructions, r_dir, float_id=None, cycles=None):
    """Plan every selected R-file/profile; bad wins, others retain existing QC.

    No rejection instruction means retain, not an upgrade to QC 1. Optional
    float/cycle filters restrict both the input instructions and output set.
    """
    from netCDF4 import Dataset
    from .writer import validate_target

    grouped = {}
    hashes = {}
    for item in instructions:
        target = item.target
        if target.selection != 'whole_profile' or target.parameters != CORE_PARAMETERS or target.sample_indices is not None or target.pressure_range is not None:
            raise ValueError('Only whole-profile core-parameter decisions are implemented')
        if item.action not in ('flag', 'no_finding') or (item.action == 'flag' and item.flag != '4'):
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
                decisions.append(Decision(target, hashes[path.name], tuple(grouped.get(target, ()))))
    return decisions
