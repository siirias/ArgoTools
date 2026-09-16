"""Resolve partner suggestions without performing scientific checks."""
from pathlib import Path

from .download import file_identity
from .instructions import CORE_PARAMETERS, Decision, sha256


def combine_instructions(instructions, r_dir):
    """Bad wins over no finding. Preserve all evidence; never infer good QC."""
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
        file_identity(target.source)
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
    decisions = []
    for target, suggestions in sorted(grouped.items(), key=lambda pair: (pair[0].source, pair[0].profile_index)):
        if any(item.action == 'flag' for item in suggestions):
            decisions.append(Decision(target, hashes[target.source], tuple(suggestions)))
    return decisions
