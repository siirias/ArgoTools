"""Shared local paths and default float for DMQC command-line tools."""
import os
from pathlib import Path, PureWindowsPath

from .instructions import read_yaml

LOCAL_SETTINGS = Path(__file__).resolve().parent.parent / 'config' / 'local.yaml'


def add_settings_argument(parser):
    parser.add_argument('--settings', type=Path, default=None,
                        help='Local settings YAML (default: config/local.yaml beside the scripts)')


def load_settings(path=None):
    """Read on demand so imports, --help, and explicit paths need no local file."""
    path = Path(path) if path is not None else LOCAL_SETTINGS
    if not path.is_file():
        raise ValueError(f'Local settings not found: {path}. Copy config/local.example.yaml '
                         'to config/local.yaml and edit it, or supply a directory explicitly.')
    document = read_yaml(path)
    if set(document) != {'work_directory', 'default_float'}:
        raise ValueError(f'{path}: expected work_directory and default_float')
    value = document['work_directory']
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{path}: work_directory must be a nonempty path string')
    # Do not silently create a directory named C:\... when running under Linux/WSL.
    if os.name != 'nt' and (PureWindowsPath(value).drive or '\\' in value):
        raise ValueError(f'{path}: use a Linux/WSL path here, e.g. /mnt/c/Data/DMQCprocessing, '
                         'rather than a Windows path')
    work_dir = Path(value).expanduser()
    if not work_dir.is_absolute():
        work_dir = path.resolve().parent / work_dir
    float_id = document['default_float']
    if type(float_id) not in (str, int) or not str(float_id).isascii() or not str(float_id).isdigit():
        raise ValueError(f'{path}: default_float must contain only digits')
    return work_dir.resolve(), str(float_id)


def processing_defaults(work_dir=None, float_id=None, settings_path=None):
    """Explicit command-line arguments take precedence over local settings."""
    if work_dir is None or float_id is None:
        configured_dir, configured_float = load_settings(settings_path)
        work_dir = configured_dir if work_dir is None else work_dir
        float_id = configured_float if float_id is None else float_id
    return Path(work_dir).expanduser(), str(float_id)


def float_directory(directory=None, settings_path=None):
    """An explicit float/R/D directory bypasses the local defaults entirely."""
    if directory is not None:
        return Path(directory).expanduser()
    work_dir, float_id = processing_defaults(settings_path=settings_path)
    return work_dir / float_id
