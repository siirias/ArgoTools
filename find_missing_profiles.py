# -*- coding: utf-8 -*-
"""
find_missing_profiles.py

Compare two .mat files (as produced by save_to_mat() in your ICES
processing script) and find profiles present in the OLD file that have
no matching counterpart in the NEW file.

A profile in OLD is considered "matched" in NEW if there's a NEW profile
within:
    - MAX_DIST_NM  nautical miles, AND
    - MAX_TIME_HOURS hours

Unmatched OLD profiles are plotted (on a coastline map if cartopy is
available) and some summary statistics are printed.

Usage:
    python find_missing_profiles.py old_file.mat new_file.mat

    or edit OLD_FILE / NEW_FILE below and run directly.
"""

import sys
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt

# ---- Edit these if running without command-line args ----------------------
OLD_FILE = r'C:\Data\DMQC\UPDATE_test\23_fmi_ctd_1501.mat'
NEW_FILE = r'C:\Data\DMQC\UPDATE_test\26v2_fmi_ctd_1501.mat'

MAX_DIST_NM = 0.5      # nautical miles
MAX_TIME_HOURS = 6.0    # hours

# Set to a longitude value (e.g. 12.0) to only include points EAST of it
# (i.e. clip off the erroneous westward extent). Set to None to disable.
WEST_LON_LIMIT = 12.0
# -----------------------------------------------------------------------------

EARTH_RADIUS_KM = 6371.0088
NM_TO_KM = 1.852


def load_profiles(mat_path):
    """Load LONG/LAT/DATES (and SOURCE if present) from a save_to_mat() file."""
    mat = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    for key in ('LONG', 'LAT', 'DATES'):
        if key not in mat:
            raise KeyError(
                f"'{key}' not found in {mat_path}. "
                f"Available keys: {[k for k in mat.keys() if not k.startswith('__')]}"
            )
    lon = np.atleast_1d(np.array(mat['LONG'], dtype=float).squeeze())
    lat = np.atleast_1d(np.array(mat['LAT'], dtype=float).squeeze())
    dates_raw = np.atleast_1d(np.array(mat['DATES']).squeeze())

    # DATES was saved as float YYYYMMDDHHMMSS -> parse to real timestamps
    dates_str = pd.Series(dates_raw).astype(float).astype('int64').astype(str)
    dates = pd.to_datetime(dates_str, format='%Y%m%d%H%M%S', errors='coerce')

    source = None
    if 'SOURCE' in mat:
        source = np.atleast_1d(np.array(mat['SOURCE']).squeeze()).astype(str)

    df = pd.DataFrame({'lon': lon, 'lat': lat, 'time': dates})
    if source is not None:
        df['source'] = source
    return df


def clip_west(df, west_lon_limit):
    """Keep only points with lon >= west_lon_limit. No-op if limit is None."""
    if west_lon_limit is None:
        return df
    before = len(df)
    clipped = df[df['lon'] >= west_lon_limit].reset_index(drop=True)
    print(f"  clipped at lon >= {west_lon_limit}: {before} -> {len(clipped)} points")
    return clipped


def find_missing(old_df, new_df, max_dist_nm=MAX_DIST_NM, max_time_hours=MAX_TIME_HOURS):
    """Return boolean mask over old_df: True where the old profile has NO
    matching profile in new_df within max_dist_nm and max_time_hours."""

    max_dist_km = max_dist_nm * NM_TO_KM
    radius_rad = max_dist_km / EARTH_RADIUS_KM
    max_time_delta = pd.Timedelta(hours=max_time_hours)

    old_lat_rad = np.radians(old_df['lat'].values)
    old_lon_rad = np.radians(old_df['lon'].values)
    new_lat_rad = np.radians(new_df['lat'].values)
    new_lon_rad = np.radians(new_df['lon'].values)

    old_time = old_df['time'].values
    new_time = new_df['time'].values

    n_old = len(old_df)
    matched = np.zeros(n_old, dtype=bool)

    try:
        from sklearn.neighbors import BallTree

        new_coords = np.column_stack([new_lat_rad, new_lon_rad])
        old_coords = np.column_stack([old_lat_rad, old_lon_rad])
        tree = BallTree(new_coords, metric='haversine')
        candidates = tree.query_radius(old_coords, r=radius_rad)

        for i, cand_idx in enumerate(candidates):
            if len(cand_idx) == 0:
                continue
            time_diffs = np.abs(new_time[cand_idx] - old_time[i])
            matched[i] = np.any(time_diffs <= max_time_delta)

    except ImportError:
        print("NOTE: scikit-learn not installed -- using slower chunked "
              "fallback. Install with: pip install scikit-learn")
        chunk_size = 500
        for start in range(0, n_old, chunk_size):
            end = min(start + chunk_size, n_old)
            lat1 = old_lat_rad[start:end][:, None]
            lon1 = old_lon_rad[start:end][:, None]
            lat2 = new_lat_rad[None, :]
            lon2 = new_lon_rad[None, :]

            dlat = lat2 - lat1
            dlon = lon2 - lon1
            a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
            dist_rad = 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

            within_dist = dist_rad <= radius_rad

            t1 = old_time[start:end][:, None]
            t2 = new_time[None, :]
            time_diff = np.abs(t2 - t1)
            within_time = time_diff <= max_time_delta

            matched[start:end] = np.any(within_dist & within_time, axis=1)

    return ~matched  # True = missing (no match found)


def print_stats(old_df, missing_mask):
    n_old = len(old_df)
    n_missing = missing_mask.sum()
    n_matched = n_old - n_missing

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Old profiles total:     {n_old}")
    print(f"Matched in new:         {n_matched} ({100 * n_matched / n_old:.1f}%)")
    print(f"Missing from new:       {n_missing} ({100 * n_missing / n_old:.1f}%)")

    missing_df = old_df[missing_mask]
    if len(missing_df) == 0:
        print("\nNo missing profiles -- everything in OLD has a match in NEW.")
        return missing_df

    # Breakdown by year
    print("\nMissing profiles by year:")
    year_counts = missing_df['time'].dt.year.value_counts().sort_index()
    for year, count in year_counts.items():
        print(f"  {year}: {count}")

    # Breakdown by rough region (1-degree lat/lon bins), top 10
    print("\nMissing profiles by 1x1 degree bin (top 10):")
    binned = missing_df.copy()
    binned['lat_bin'] = np.floor(binned['lat'])
    binned['lon_bin'] = np.floor(binned['lon'])
    bin_counts = binned.groupby(['lat_bin', 'lon_bin']).size().sort_values(ascending=False)
    for (lat_bin, lon_bin), count in bin_counts.head(10).items():
        print(f"  lat {lat_bin:.0f}-{lat_bin+1:.0f}, lon {lon_bin:.0f}-{lon_bin+1:.0f}: {count}")

    # Bounding box of missing points, for reference
    print(f"\nMissing points bounding box:")
    print(f"  lon: {missing_df['lon'].min():.2f} to {missing_df['lon'].max():.2f}")
    print(f"  lat: {missing_df['lat'].min():.2f} to {missing_df['lat'].max():.2f}")

    if 'source' in missing_df.columns:
        print("\nSample of missing SOURCE ids (up to 10):")
        for s in missing_df['source'].head(10):
            print(f"  {s}")

    return missing_df


def diagnose_near_misses(old_df, new_df, missing_mask,
                          max_dist_nm=MAX_DIST_NM, max_time_hours=MAX_TIME_HOURS):
    """For each MISSING old profile, find its single nearest new profile in
    space (ignoring the time window entirely) and report how close it
    actually is, plus the time gap to that nearest neighbor. This tells you
    whether "missing" mostly means:
      - a near-identical station/time exists but just outside the time
        window  -> timestamp handling issue
      - a near-identical station exists at almost the same spot but outside
        the distance threshold -> coordinate precision/rounding issue
      - nothing close exists at all -> genuinely absent from the new file
    """
    missing_df = old_df[missing_mask]
    if len(missing_df) == 0:
        return

    lat1 = np.radians(missing_df['lat'].values)
    lon1 = np.radians(missing_df['lon'].values)
    lat2 = np.radians(new_df['lat'].values)
    lon2 = np.radians(new_df['lon'].values)

    try:
        from sklearn.neighbors import BallTree
        new_coords = np.column_stack([lat2, lon2])
        old_coords = np.column_stack([lat1, lon1])
        tree = BallTree(new_coords, metric='haversine')
        dist_rad, idx = tree.query(old_coords, k=1)
        dist_km = dist_rad[:, 0] * EARTH_RADIUS_KM
        nearest_idx = idx[:, 0]
    except ImportError:
        # Chunked fallback
        n = len(missing_df)
        dist_km = np.full(n, np.nan)
        nearest_idx = np.full(n, -1, dtype=int)
        chunk_size = 500
        for start in range(0, n, chunk_size):
            end = min(start + chunk_size, n)
            la1 = lat1[start:end][:, None]
            lo1 = lon1[start:end][:, None]
            dlat = lat2[None, :] - la1
            dlon = lon2[None, :] - lo1
            a = np.sin(dlat / 2) ** 2 + np.cos(la1) * np.cos(lat2[None, :]) * np.sin(dlon / 2) ** 2
            d = 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))) * EARTH_RADIUS_KM
            best = np.argmin(d, axis=1)
            dist_km[start:end] = d[np.arange(end - start), best]
            nearest_idx[start:end] = best

    nm_thresh_km = max_dist_nm * NM_TO_KM
    new_time = new_df['time'].values
    old_time = missing_df['time'].values
    time_diff_hours = np.abs((new_time[nearest_idx] - old_time)) / np.timedelta64(1, 'h')

    print("\n" + "=" * 60)
    print(f"NEAR-MISS DIAGNOSTIC (nearest new profile, any time)")
    print("=" * 60)

    bins_km = [0, 0.1, nm_thresh_km, 1, 5, 20, np.inf]
    labels = ['<0.1km (~same spot)', f'0.1km-{nm_thresh_km:.2f}km (within nm thresh)',
              f'{nm_thresh_km:.2f}-1km', '1-5km', '5-20km', '>20km']
    cats = pd.cut(dist_km, bins=bins_km, labels=labels, right=False)
    print("\nDistance to nearest new profile (ignoring time entirely):")
    for label, count in cats.value_counts().reindex(labels).items():
        pct = 100 * count / len(missing_df)
        print(f"  {label:35s}: {count:6d} ({pct:5.1f}%)")

    # Of those that ARE within the distance threshold spatially,
    # how far off in time are they? Use a distribution, not just a median --
    # for recurring/fixed monitoring stations, "nearest in space" often finds
    # a visit from a totally different year, which would badly skew a median
    # and wrongly suggest these are near-misses on the time window.
    within_dist = dist_km <= nm_thresh_km
    n_within_dist = within_dist.sum()
    if n_within_dist > 0:
        print(f"\nOf the {n_within_dist} missing profiles with a spatially close "
              f"new profile (<{max_dist_nm} nm) but still unmatched,")
        print(f"time gap to that nearest-in-space profile, bucketed:")
        td = time_diff_hours[within_dist]
        time_bins_h = [0, max_time_hours, 24, 24 * 30, 24 * 365, np.inf]
        time_labels = [f'<{max_time_hours:.0f}h (should have matched!)',
                       f'{max_time_hours:.0f}h-1day', '1day-30days',
                       '30days-1yr', '>1yr (different visit, same station)']
        td_cats = pd.cut(td, bins=time_bins_h, labels=time_labels, right=False)
        for label, count in td_cats.value_counts().reindex(time_labels).items():
            pct = 100 * count / n_within_dist
            print(f"  {label:38s}: {count:6d} ({pct:5.1f}%)")
        n_true_near_miss = (td <= 24 * 30).sum()
        print(f"\n  -> {n_true_near_miss} of these ({100*n_true_near_miss/n_within_dist:.1f}%) are true "
              f"near-misses (within ~1 month) that a looser time window could recover.")
        print(f"  -> the rest are the SAME station being visited at a totally different time -- "
              f"the station exists in the new file, just not this specific date's profile.")

    n_far = (dist_km > 20).sum()
    print(f"\n{n_far} missing profiles ({100*n_far/len(missing_df):.1f}%) have no new "
          f"profile within 20 km at all -> likely genuinely absent from the new file.")


def diagnose_by_cruise(old_df, missing_mask, min_count=5, top_n=20):
    """Group profiles by cruise id (parsed from SOURCE, e.g. 'ices_067I_0079'
    -> cruise '067I') and report which cruises are mostly/entirely missing
    from the new file. A whole cruise missing (rather than scattered single
    stations) points to an entire campaign being dropped, not a matching
    threshold issue."""
    if 'source' not in old_df.columns:
        print("\n(No SOURCE column available -- skipping per-cruise diagnostic.)")
        return

    df = old_df.copy()
    df['missing'] = missing_mask
    # SOURCE looks like 'ices_<cruise>_<station>' -> take the cruise part
    df['cruise'] = df['source'].str.split('_').str[1]
    df['year'] = df['time'].dt.year

    grouped = df.groupby('cruise').agg(
        sum=('missing', 'sum'),
        count=('missing', 'count'),
        year_min=('year', 'min'),
        year_max=('year', 'max'),
    )
    grouped['pct_missing'] = 100 * grouped['sum'] / grouped['count']
    grouped = grouped[grouped['count'] >= min_count]
    grouped = grouped.sort_values(['pct_missing', 'count'], ascending=[False, False])

    def year_label(row):
        if row['year_min'] == row['year_max']:
            return f"{int(row['year_min'])}"
        return f"{int(row['year_min'])}-{int(row['year_max'])}"

    print("\n" + "=" * 60)
    print(f"BY-CRUISE DIAGNOSTIC (cruises with >= {min_count} old profiles)")
    print("=" * 60)
    print(f"{'cruise':10s} {'year(s)':12s} {'missing':>8s} {'total':>8s} {'% missing':>10s}")
    for cruise, row in grouped.head(top_n).iterrows():
        print(f"{str(cruise):10s} {year_label(row):12s} {int(row['sum']):8d} "
              f"{int(row['count']):8d} {row['pct_missing']:9.1f}%")

    fully_missing = grouped[grouped['pct_missing'] >= 99.9]
    if len(fully_missing) > 0:
        print(f"\n{len(fully_missing)} cruise(s) are essentially 100% missing "
              f"({fully_missing['sum'].sum()} profiles total):")
        for cruise, row in fully_missing.iterrows():
            print(f"  {cruise} ({year_label(row)}): {int(row['sum'])} profiles")
        print("  -> these look like entire campaigns dropped from the new export, "
              "not a matching/threshold artifact.")


def plot_missing(old_df, new_df, missing_mask):
    missing_df = old_df[missing_mask]

    pad = 0.5
    all_lon = pd.concat([old_df['lon'], new_df['lon']])
    all_lat = pd.concat([old_df['lat'], new_df['lat']])
    lon_min, lon_max = all_lon.min() - pad, all_lon.max() + pad
    lat_min, lat_max = all_lat.min() - pad, all_lat.max() + pad

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

        # faint context: all new points
        ax.scatter(new_df['lon'], new_df['lat'], s=8, c='lightgray', alpha=0.4,
                   label=f'New (all, {len(new_df)} pts)', transform=ccrs.PlateCarree(), zorder=1)
        ax.scatter(missing_df['lon'], missing_df['lat'], s=25, c='red', alpha=0.8,
                   label=f'Missing from new ({len(missing_df)} pts)',
                   transform=ccrs.PlateCarree(), zorder=3)

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

        ax.scatter(new_df['lon'], new_df['lat'], s=8, c='lightgray', alpha=0.4,
                   label=f'New (all, {len(new_df)} pts)', zorder=1)
        ax.scatter(missing_df['lon'], missing_df['lat'], s=25, c='red', alpha=0.8,
                   label=f'Missing from new ({len(missing_df)} pts)', zorder=3)

    ax.legend(loc='upper right')
    ax.set_title(f'Old profiles with no match in new\n'
                 f'(within {MAX_DIST_NM} nm & {MAX_TIME_HOURS} h)')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) == 3:
        old_file, new_file = sys.argv[1], sys.argv[2]
    else:
        old_file, new_file = OLD_FILE, NEW_FILE
        print(f"No command-line args given, using defaults:\n  {old_file}\n  {new_file}")

    old_df = load_profiles(old_file)
    new_df = load_profiles(new_file)

    print(f"Loaded {len(old_df)} old profiles, {len(new_df)} new profiles.")

    if WEST_LON_LIMIT is not None:
        print(f"\nApplying west clip at lon >= {WEST_LON_LIMIT}:")
        old_df = clip_west(old_df, WEST_LON_LIMIT)
        new_df = clip_west(new_df, WEST_LON_LIMIT)

    print(f"\nMatching within {MAX_DIST_NM} nm and {MAX_TIME_HOURS} hours...\n")

    missing_mask = find_missing(old_df, new_df)
    missing_df = print_stats(old_df, missing_mask)
    diagnose_near_misses(old_df, new_df, missing_mask)
    diagnose_by_cruise(old_df, missing_mask)
    plot_missing(old_df, new_df, missing_mask)