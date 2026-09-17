"""Inspect source uncertainties and write constant-per-parameter instructions."""
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from .download import file_identity
from .instructions import CORE_PARAMETERS, sha256, uncertainty_value, write_instructions

UNITS = {'PRES': 'decibar', 'TEMP': 'degree_Celsius', 'PSAL': 'psu'}


def inspect_sources(directory):
    """Read a float folder or R/ folder; resolution is never an uncertainty default."""
    directory = Path(directory)
    r_dir = directory / 'R' if (directory / 'R').is_dir() else directory
    paths = sorted(r_dir.glob('R*.nc'))
    if not paths:
        raise ValueError(f'No R*.nc files in {r_dir}')
    floats = {file_identity(path.name)[0] for path in paths}
    if len(floats) != 1:
        raise ValueError('Select R-files from exactly one float')
    stats = {p: {'count': 0, 'minimum': None, 'maximum': None, 'units': UNITS[p]} for p in CORE_PARAMETERS}
    sources = []
    for path in paths:
        checksum = sha256(path)
        with Dataset(path) as ds:
            count = len(ds.dimensions['N_PROF'])
            for p in CORE_PARAMETERS:
                var = ds[p + '_ADJUSTED_ERROR']
                units = getattr(var, 'units', None)
                if units != UNITS[p]:
                    raise ValueError(f'{path.name}: {var.name} units must be {UNITS[p]!r}, got {units!r}')
                values = np.ma.masked_invalid(var[:])
                flags = np.ma.filled(ds[p + '_ADJUSTED_QC'][:], b' ')
                values = np.ma.masked_where(~np.isin(flags, [b'1', b'2', b'3', b'5', b'8']), values).compressed()
                values = values[values > 0]
                if values.size:
                    stat = stats[p]
                    lo, hi = float(values.min()), float(values.max())
                    stat['minimum'] = lo if stat['minimum'] is None else min(stat['minimum'], lo)
                    stat['maximum'] = hi if stat['maximum'] is None else max(stat['maximum'], hi)
                    stat['count'] += int(values.size)
        if checksum != sha256(path):
            raise ValueError(f'Source changed while reading: {path}')
        sources.append({'source': path.name, 'profiles': count, 'source_sha256': checksum})
    for stat in stats.values():
        # Float32 storage adds insignificant digits; display a useful suggested value.
        stat['default'] = (float(format(stat['minimum'], '.7g'))
                           if stat['minimum'] is not None and stat['minimum'] == stat['maximum'] else None)
    return r_dir, sources, stats


def write_defaults(path, sources, values, reasons=None):
    """Write float-wide defaults, grouping parameters with the same justification."""
    if any(p not in CORE_PARAMETERS for p in values):
        raise ValueError('Only PRES, TEMP and PSAL are supported')
    values = {p: uncertainty_value(v) for p, v in values.items()}
    floats = {file_identity(source['source'])[0] for source in sources}
    if len(floats) != 1:
        raise ValueError('Defaults require sources from exactly one float')
    float_id = floats.pop()
    groups = {}
    for p, value in values.items():
        reason = (reasons or {}).get(p, 'Operator-selected default')
        groups.setdefault(reason, {})[p] = value
    entries = [{'target': {'float': float_id, 'selection': 'all_profiles'},
                'action': 'set_uncertainty', 'values': grouped, 'reason': reason}
               for reason, grouped in groups.items()]
    return write_instructions(path, 'default_uncertainties', entries)
