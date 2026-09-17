"""Read the parameter inventory of an individual Argo N_PROF entry."""
import numpy as np

from .instructions import CORE_PARAMETERS


def profile_parameters(ds, index):
    """Return declared parameters, rejecting malformed or ambiguous inventories.

    Empty padded slots are ignored. Presence means declared in STATION_PARAMETERS,
    not that a parameter happens to have a usable value in this particular profile.
    """
    if 'STATION_PARAMETERS' not in ds.variables:
        raise ValueError('Missing STATION_PARAMETERS')
    var = ds['STATION_PARAMETERS']
    if var.dimensions != ('N_PROF', 'N_PARAM', 'STRING16'):
        raise ValueError('Unsupported STATION_PARAMETERS dimensions')
    # Read S1 explicitly even if a caller enables netCDF string conversion.
    var.set_auto_chartostring(False)
    names = tuple(np.ma.filled(row, b' ').tobytes().decode('ascii').strip(' \x00')
                  for row in var[index])
    names = tuple(name for name in names if name)
    if not names or len(set(names)) != len(names):
        raise ValueError(f'Profile {index}: STATION_PARAMETERS must be nonempty and unique')
    return names


def core_parameters(ds, index):
    names = profile_parameters(ds, index)
    return tuple(p for p in CORE_PARAMETERS if p in names)
