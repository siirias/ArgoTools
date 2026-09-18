"""Fetch immutable local R-file snapshots; no QC decisions are made here."""
import os
import re
import tempfile
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlopen


FILE_PATTERN = re.compile(r"R(?P<float>\d+)_(?P<cycle>\d+)(?P<direction>D?)\.nc")


def file_identity(name):
    match = FILE_PATTERN.fullmatch(name)
    if not match:
        raise ValueError(f"Expected an Argo R-file name, got {name!r}")
    return match['float'], int(match['cycle'])


def download_r_files(float_id, r_dir, dac='coriolis', cycles=None, dry_run=False,
                     base_root='https://data-argo.ifremer.fr/dac'):
    """Skip existing sources; publish a download only after it is complete."""
    from netCDF4 import Dataset

    if not str(float_id).isdigit() or not re.fullmatch(r'[A-Za-z0-9_-]+', dac):
        raise ValueError('Invalid float ID or DAC')
    directory = Path(r_dir)
    base_url = f'{base_root}/{dac}/{float_id}/profiles/'
    with urlopen(base_url, timeout=60) as response:
        html = response.read().decode('utf-8', errors='replace')
    names = sorted(set(re.findall(r'href=[\"\']([^\"\']+)[\"\']', html)))
    downloaded = []
    for name in names:
        match = FILE_PATTERN.fullmatch(name)
        if not match or match['float'] != str(float_id):
            continue
        if cycles is not None and int(match['cycle']) not in cycles:
            continue
        path = directory / name
        if path.exists():
            continue
        print(f'{"Would download" if dry_run else "Downloading"}: {name}')
        if dry_run:
            continue
        directory.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(dir=directory, suffix='.part')
        try:
            with os.fdopen(fd, 'wb') as output, urlopen(urljoin(base_url, name), timeout=60) as response:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            with Dataset(temporary):
                pass
            # Hard-link publication refuses to replace an existing source.
            os.link(temporary, path)
            downloaded.append(path)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return downloaded
