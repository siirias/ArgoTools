# -*- coding: utf-8 -*-
"""
compare_mat_points.py

Plot LONG/LAT points from two .mat files (as produced by save_to_mat() in
your ICES processing script) on one map, so you can visually compare
coverage between an "old" and a "new" dataset.

Usage:
    python compare_mat_points.py old_file.mat new_file.mat

    or just edit OLD_FILE / NEW_FILE below and run directly.
"""

import sys
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

# ---- Edit these if running without command-line args ----------------------
OLD_FILE = r'C:\Data\DMQC\UPDATE_test\23_fmi_ctd_1501.mat'
NEW_FILE = r'C:\Data\DMQC\UPDATE_test\26v2_fmi_ctd_1501.mat'
# -----------------------------------------------------------------------------

def load_lonlat(mat_path):
    """Load LONG/LAT arrays from a .mat file saved by save_to_mat()."""
    mat = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    if 'LONG' not in mat or 'LAT' not in mat:
        raise KeyError(
            f"'LONG'/'LAT' not found in {mat_path}. "
            f"Available keys: {[k for k in mat.keys() if not k.startswith('__')]}"
        )
    lon = np.atleast_1d(np.array(mat['LONG'], dtype=float).squeeze())
    lat = np.atleast_1d(np.array(mat['LAT'], dtype=float).squeeze())
    return lon, lat


def plot_comparison(old_file, new_file):
    lon_old, lat_old = load_lonlat(old_file)
    lon_new, lat_new = load_lonlat(new_file)

    print(f"Old file: {old_file}  -> {len(lon_old)} points")
    print(f"New file: {new_file} -> {len(lon_new)} points")

    # Bounding box with a little padding, based on both datasets combined
    all_lon = np.concatenate([lon_old, lon_new])
    all_lat = np.concatenate([lat_old, lat_new])
    pad = 0.5
    lon_min, lon_max = np.nanmin(all_lon) - pad, np.nanmax(all_lon) + pad
    lat_min, lat_max = np.nanmin(all_lat) - pad, np.nanmax(all_lat) + pad

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        fig = plt.figure(figsize=(10, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

        ax.add_feature(cfeature.LAND, facecolor='#e8e4d8', zorder=0)
        ax.add_feature(cfeature.OCEAN, facecolor='#dbeaf2', zorder=0)
        ax.coastlines(resolution='10m', linewidth=0.8, zorder=1)
        ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

        ax.scatter(lon_old, lat_old, s=25, c='blue', alpha=0.2,
                   label=f'Old ({len(lon_old)} pts)', transform=ccrs.PlateCarree(),
                   zorder=2)
        ax.scatter(lon_new, lat_new, s=15, c='red', alpha=0.2,
                   label=f'New ({len(lon_new)} pts)', transform=ccrs.PlateCarree(),
                   zorder=3)

    except ImportError:
        print("NOTE: cartopy not installed -- plotting plain lon/lat scatter "
              "without coastlines. Install with: pip install cartopy")
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_aspect('equal')
        ax.grid(True, linewidth=0.3, alpha=0.5)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        ax.scatter(lon_old, lat_old, s=25, c='blue', alpha=0.6,
                   label=f'Old ({len(lon_old)} pts)', zorder=2)
        ax.scatter(lon_new, lat_new, s=15, c='red', alpha=0.8,
                   label=f'New ({len(lon_new)} pts)', zorder=3)

    ax.legend(loc='upper right')
    ax.set_title('Point coverage comparison: Old (blue) vs New (red)')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) == 3:
        old_file, new_file = sys.argv[1], sys.argv[2]
    else:
        old_file, new_file = OLD_FILE, NEW_FILE
        print(f"No command-line args given, using defaults:\n  {old_file}\n  {new_file}")

    plot_comparison(old_file, new_file)