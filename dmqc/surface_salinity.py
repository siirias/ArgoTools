"""Regional surface salinity check; raw data in, whole-profile instructions out."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import math

import numpy as np
from netCDF4 import Dataset

from .download import file_identity
from .instructions import read_yaml, sha256, write_instructions


def number(value):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f'Expected a finite number, got {value!r}')
    return float(value)


def keys(value, allowed, required):
    if not isinstance(value, dict) or set(value) - set(allowed) or set(required) - set(value):
        raise ValueError(f'Expected mapping with required keys {required}; allowed keys {allowed}')


def interval(value, limits):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError('Ranges must be [minimum, maximum]')
    lo, hi = map(number, value)
    if not limits[0] <= lo < hi <= limits[1]:
        raise ValueError(f'Invalid range {value} within {limits}')
    return lo, hi


def load_config(path):
    config = read_yaml(path)
    keys(config, ('surface', 'areas'), ('surface', 'areas'))
    surface = config['surface']
    fields = ('pressure_range_dbar', 'statistic', 'minimum_samples', 'pressure_qc', 'salinity_qc')
    keys(surface, fields, fields)
    interval(surface['pressure_range_dbar'], (-100, 12000))
    if surface['statistic'] not in ('median', 'maximum'):
        raise ValueError('statistic must be median or maximum')
    if type(surface['minimum_samples']) is not int or surface['minimum_samples'] < 1:
        raise ValueError('minimum_samples must be a positive integer')
    for field in ('pressure_qc', 'salinity_qc'):
        flags = surface[field]
        if not isinstance(flags, list) or not flags or any(type(q) is not str or q not in tuple('01234589') for q in flags):
            raise ValueError(f'{field}: expected quoted QC codes, e.g. ["1", "4"]')
    if not isinstance(config['areas'], list) or not config['areas']:
        raise ValueError('areas must be a nonempty list')
    names = set()
    for area in config['areas']:
        fields = ('name', 'latitude', 'longitude', 'maximum_surface_salinity')
        keys(area, fields, fields)
        name = area['name']
        if not isinstance(name, str) or not name.strip() or name in names:
            raise ValueError('Area names must be nonempty and unique')
        names.add(name)
        if number(area['maximum_surface_salinity']) < 0:
            raise ValueError('Salinity limits must be nonnegative')
        interval(area['latitude'], (-90, 90))
        interval(area['longitude'], (-180, 180))
    return config


def matching_areas(config, latitude, longitude):
    """Rectangles include their lower edges and exclude their upper edges."""
    return [area for area in config['areas']
            if area['latitude'][0] <= latitude < area['latitude'][1]
            and area['longitude'][0] <= longitude < area['longitude'][1]]


def coordinate(ds, name, ip):
    value = ds[name][ip]
    if np.ma.is_masked(value) or not np.isfinite(value):
        return None
    return float(value)


def check_directory(directory, config):
    """Read every source profile independently. A malformed file aborts saving."""
    root = Path(directory)
    r_dir = root / 'R' if (root / 'R').is_dir() else root
    paths = sorted(r_dir.glob('R*.nc'))
    if not paths:
        raise ValueError(f'No R*.nc files in {r_dir}')
    if len({file_identity(p.name)[0] for p in paths}) != 1:
        raise ValueError('Select files from one float')
    surface = config['surface']
    results, instructions, sources = [], [], []
    for path in paths:
        checksum = sha256(path)
        with Dataset(path) as ds:
            ds.set_auto_chartostring(False)
            for name, dims in {'LATITUDE': ('N_PROF',), 'LONGITUDE': ('N_PROF',),
                               'PRES': ('N_PROF', 'N_LEVELS'), 'PSAL': ('N_PROF', 'N_LEVELS'),
                               'PRES_QC': ('N_PROF', 'N_LEVELS'), 'PSAL_QC': ('N_PROF', 'N_LEVELS')}.items():
                if name not in ds.variables or ds[name].dimensions != dims:
                    raise ValueError(f'{path.name}: missing or unsupported {name}')
            for ip in range(len(ds.dimensions['N_PROF'])):
                lat, lon = coordinate(ds, 'LATITUDE', ip), coordinate(ds, 'LONGITUDE', ip)
                result = {'source': path.name, 'profile_index': ip, 'latitude': lat,
                          'longitude': lon, 'status': 'not_evaluated'}
                results.append(result)
                if lat is None or lon is None or not -90 <= lat <= 90 or not -180 <= lon <= 180:
                    result['reason'] = 'missing_or_invalid_position'
                    continue
                areas = matching_areas(config, lat, lon)
                if not areas:
                    result['reason'] = 'no_matching_area'
                    continue
                pressure = np.ma.masked_invalid(ds['PRES'][ip])
                salinity = np.ma.masked_invalid(ds['PSAL'][ip])
                pq = np.ma.filled(ds['PRES_QC'][ip], b' ')
                sq = np.ma.filled(ds['PSAL_QC'][ip], b' ')
                lo, hi = surface['pressure_range_dbar']
                eligible = (~np.ma.getmaskarray(pressure) & ~np.ma.getmaskarray(salinity)
                            & np.ma.filled((pressure >= lo) & (pressure <= hi), False)
                            & np.isin(pq, [q.encode() for q in surface['pressure_qc']])
                            & np.isin(sq, [q.encode() for q in surface['salinity_qc']]))
                values = salinity.data[eligible]
                result['sample_count'] = int(values.size)
                result['matching_areas'] = [a['name'] for a in areas]
                if values.size < surface['minimum_samples']:
                    result['reason'] = 'insufficient_surface_samples'
                    continue
                value = float(np.median(values) if surface['statistic'] == 'median' else np.max(values))
                result['surface_salinity'] = value
                result['checks'] = [{'area': a['name'], 'limit': a['maximum_surface_salinity'],
                                     'exceeded': value > a['maximum_surface_salinity']} for a in areas]
                failed = [a for a in result['checks'] if a['exceeded']]
                result['status'] = 'reject' if failed else 'pass'
                if failed:
                    reason = '; '.join(f"{a['area']}: surface {surface['statistic']} PSAL {value:g} > {a['limit']:g}" for a in failed)
                    reason += f' ({lo:g}–{hi:g} dbar, {values.size} samples)'
                    instructions.append({'target': {'source': path.name, 'profile_index': ip,
                                                     'selection': 'whole_profile'},
                                         'action': 'flag', 'flag': '4', 'reason': reason,
                                         'source_sha256': checksum})
        if checksum != sha256(path):
            raise ValueError(f'Source changed during checking: {path}')
        sources.append({'source': path.name, 'sha256': checksum})
    report = {'schema_version': 1, 'checker': 'surface_salinity',
              'created_utc': datetime.now(timezone.utc).isoformat(), 'configuration': config,
              'sources': sources, 'summary': dict(Counter(r['status'] for r in results)),
              'skip_reasons': dict(Counter(r['reason'] for r in results if r['status'] == 'not_evaluated')),
              'profiles': results}
    return r_dir, instructions, report


def save_instructions(path, entries):
    return write_instructions(path, 'surface_salinity', entries)
