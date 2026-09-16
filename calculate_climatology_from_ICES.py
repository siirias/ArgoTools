import argparse

import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path
import matplotlib.pyplot as plt

POINT_COLUMNS = ["LONG", "LAT", "PRES", "value"]
STAT_COLUMNS = ["PRES_window", "count", "mean", "std", "min", "25%", "50%", "75%", "max", "midbin"]
# Historical Bothnian Sea QC rule: reject all salinity > 20, including the
# known placeholder value 35. This is an area-specific upper limit.
BOTHSEA_MAX_SALINITY = 20


def empty_points():
    return pd.DataFrame({column: pd.Series(dtype=float) for column in POINT_COLUMNS})


def profiles_to_points(mat, value_name):
    """Flatten MATLAB cell-array profiles into finite observations.

    Input arrays must retain their MATLAB dimensions (loadmat without squeeze).
    Each cell in PRES and the requested variable contains one profile.
    """
    lon = np.asarray(mat["LONG"], dtype=float).ravel()
    lat = np.asarray(mat["LAT"], dtype=float).ravel()
    pressure = np.asarray(mat["PRES"], dtype=object).ravel()
    values = np.asarray(mat[value_name], dtype=object).ravel()
    if not (len(lon) == len(lat) == len(pressure) == len(values)):
        raise ValueError(f"Mismatched profile counts in LONG, LAT, PRES and {value_name}")

    columns = {column: [] for column in POINT_COLUMNS}
    for i, (x, y, pres, val) in enumerate(zip(lon, lat, pressure, values)):
        pres = np.asarray(pres, dtype=float).ravel()
        val = np.asarray(val, dtype=float).ravel()
        if pres.size != val.size:
            raise ValueError(f"Profile {i}: PRES and {value_name} have different lengths")
        ok = np.isfinite(pres) & np.isfinite(val) & np.isfinite(x) & np.isfinite(y)
        if ok.any():
            columns["LONG"].append(np.full(ok.sum(), x))
            columns["LAT"].append(np.full(ok.sum(), y))
            columns["PRES"].append(pres[ok])
            columns["value"].append(val[ok])

    if not columns["PRES"]:
        return empty_points()
    return pd.DataFrame({key: np.concatenate(parts) for key, parts in columns.items()})


def validate_year_range(first_year, last_year):
    for name, year in (("first_year", first_year), ("last_year", last_year)):
        if year is not None and (
            isinstance(year, bool) or not isinstance(year, (int, np.integer))
            or not 1 <= year <= 9999
        ):
            raise ValueError(f"{name} must be an integer between 1 and 9999")
    if first_year is not None and last_year is not None and first_year > last_year:
        raise ValueError("first_year must be less than or equal to last_year")


def load_filtered_profiles(matfile, value_names, first_year=None, last_year=None):
    """Filter profiles immediately after loading, using inclusive year bounds.

    ICES DATES values are numeric YYYYMMDDhhmmss timestamps. Profiles with
    nonfinite dates are excluded when either year bound is supplied.
    """
    validate_year_range(first_year, last_year)
    columns = list(dict.fromkeys(["LONG", "LAT", "PRES", *value_names]))
    bounded = first_year is not None or last_year is not None
    mat = scipy.io.loadmat(matfile, variable_names=columns + (["DATES"] if bounded else []))
    if not bounded:
        return mat
    if "DATES" not in mat:
        raise ValueError(f"{matfile}: DATES is required for year filtering")
    dates = np.asarray(mat["DATES"], dtype=float).ravel()
    years = np.floor(dates / 1e10)
    keep = np.isfinite(dates)
    if first_year is not None:
        keep &= years >= first_year
    if last_year is not None:
        keep &= years <= last_year
    for column in columns:
        values = np.asarray(mat[column])
        if values.size != dates.size:
            raise ValueError(f"{matfile}: mismatched profile counts in DATES and {column}")
        mat[column] = values.ravel()[keep].reshape(1, -1)
    return mat


def load_profiles_from_mat(matfile, value_name, first_year=None, last_year=None):
    """Load one variable, optionally restricted to an inclusive year range."""
    mat = load_filtered_profiles(matfile, [value_name], first_year, last_year)
    return profiles_to_points(mat, value_name)


def load_all_mat_points(matfiles, value_names, first_year=None, last_year=None):
    """Read each file once; return points keyed by MATLAB variable name.

    A single variable name also remains supported and returns one DataFrame.
    """
    validate_year_range(first_year, last_year)
    single_variable = isinstance(value_names, str)
    names = [value_names] if single_variable else list(value_names)
    parts = {name: [] for name in names}
    for matfile in matfiles:
        mat = load_filtered_profiles(matfile, names, first_year, last_year)
        for name in names:
            try:
                points = profiles_to_points(mat, name)
            except (KeyError, ValueError, TypeError) as exc:
                raise ValueError(f"{matfile} / {name}: {exc}") from exc
            if not points.empty:
                parts[name].append(points)
    result = {
        name: pd.concat(frames, ignore_index=True) if frames else empty_points()
        for name, frames in parts.items()
    }
    return result[names[0]] if single_variable else result


def make_statistics(df):
    """Describe observations in [0, 1], (1, 3], (3, 5], ... dbar bins.

    Every observation has equal weight: profiles with more samples in a bin
    contribute more. Empty bins are retained. Negative pressure is invalid.
    """
    if df.empty:
        return pd.DataFrame(columns=STAT_COLUMNS)
    if not np.isfinite(df[["PRES", "value"]].to_numpy(dtype=float)).all():
        raise ValueError("Statistics require finite pressure and value observations")
    if (df["PRES"] < 0).any():
        raise ValueError("Statistics require nonnegative pressure")

    edges = np.r_[0, 1, np.arange(3, np.ceil(df["PRES"].max()) + 2, 2)]
    windows = pd.cut(df["PRES"], bins=edges, include_lowest=True)
    stats = df["value"].groupby(windows, observed=False).describe()
    stats.index.name = "PRES_window"
    # include_lowest adjusts the displayed first interval slightly below zero;
    # use physical bin edges, not that display interval, for exact midpoints.
    stats["midbin"] = (edges[:-1] + edges[1:]) / 2
    return stats.reset_index()


def plot_selected_profiles(selected, area, out_path):
    """Plot accepted profile positions and the rectangular area box."""
    if selected.empty:
        return

    # One dot per unique location; repeated casts at a station share a dot.
    pts = selected[["LONG", "LAT"]].drop_duplicates()

    min_lon = area["min_lon"]
    max_lon = area["max_lon"]
    min_lat = area["min_lat"]
    max_lat = area["max_lat"]

    fig, ax = plt.subplots(figsize=(7, 7))

    ax.scatter(pts["LONG"], pts["LAT"], s=8, alpha=0.6)

    # area rectangle
    box_lon = [min_lon, max_lon, max_lon, min_lon, min_lon]
    box_lat = [min_lat, min_lat, max_lat, max_lat, min_lat]
    ax.plot(box_lon, box_lat, linestyle="--", linewidth=1.5)

    ax.set_xlim(min_lon - 0.5, max_lon + 0.5)
    ax.set_ylim(min_lat - 0.5, max_lat + 0.5)

    ax.set_xlabel("Longitude [°E]")
    ax.set_ylabel("Latitude [°N]")
    ax.set_title(f"{area['area']} accepted unique locations, n={len(pts)}")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

def plot_profile_cloud(selected, area_name, var_name, out_path):
    """Plot all accepted observations as a pressure-value cloud."""

    if selected.empty:
        return

    fig, ax = plt.subplots(figsize=(5, 8))

    ax.scatter(
        selected["value"],
        selected["PRES"],
        s=2,
        alpha=0.05,
        edgecolors="none"
    )

    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)

    ax.set_xlabel(var_name)
    ax.set_ylabel("Pressure [dbar]")
    ax.set_title(f"{area_name}: {var_name} ({len(selected)} observations)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

def main(w_dir=Path('/mnt/c/Data/DMQC/UPDATE_test/'), first_year=None, last_year=None):
    """Write statistics and plots, optionally limited to inclusive calendar years."""
    validate_year_range(first_year, last_year)
    w_dir = Path(w_dir)
    matfiles = [w_dir / f"fmi_ctd_{region}.mat" for region in (1601, 1602, 1501, 1502)]
    areas = pd.read_csv(w_dir / "Climatology_boxes.txt")
    variables = {
        "Practical_Salinity_dmnless": "SAL",
        "Temperature_degC": "TEMP",
    }
    points_by_variable = load_all_mat_points(
        matfiles, variables.values(), first_year=first_year, last_year=last_year
    )

    for out_var_name, mat_var in variables.items():
        all_points = points_by_variable[mat_var]
        for _, area in areas.iterrows():
            area_name = area["area"]
            selected = all_points[
                all_points["LONG"].between(area["min_lon"], area["max_lon"])
                & all_points["LAT"].between(area["min_lat"], area["max_lat"])
            ]
            if area_name == "BothSea" and mat_var == "SAL":
                rejected = selected["value"] > BOTHSEA_MAX_SALINITY
                print(f"{area_name} / {out_var_name}: removed {rejected.sum()} "
                      f"observations with salinity > {BOTHSEA_MAX_SALINITY}")
                selected = selected[~rejected]
            if selected.empty:
                print(f"No data for {area_name} / {out_var_name}")
                continue

            stats = make_statistics(selected)
            outfile = w_dir / f"ICES_Statistics_{out_var_name}_by_depth_in_{area_name}.csv"
            stats.to_csv(outfile)  # Retain the legacy unnamed index column.
            mapfile = w_dir / f"Map_{out_var_name}_in_{area_name}.png"
            plot_selected_profiles(selected, area, mapfile)
            profilefile = w_dir / f"Profiles_{out_var_name}_in_{area_name}.png"
            plot_profile_cloud(selected, area_name, out_var_name, profilefile)
            print(outfile, selected.shape, stats.shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--w-dir", type=Path, default=Path('/mnt/c/Data/DMQC/UPDATE_test/'))
    parser.add_argument("--first-year", type=int, help="First year to include (inclusive)")
    parser.add_argument("--last-year", type=int, help="Last year to include (inclusive)")
    args = parser.parse_args()
    try:
        validate_year_range(args.first_year, args.last_year)
    except ValueError as exc:
        parser.error(str(exc))
    main(args.w_dir, first_year=args.first_year, last_year=args.last_year)
