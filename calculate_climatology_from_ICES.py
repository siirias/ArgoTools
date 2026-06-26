import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path
import matplotlib.pyplot as plt

def load_profiles_from_mat(matfile, value_name):
    mat = scipy.io.loadmat(matfile, squeeze_me=True)

    rows = []
    for lon, lat, pres, val in zip(mat["LONG"], mat["LAT"], mat["PRES"], mat[value_name]):
        pres = np.asarray(pres, dtype=float).ravel()
        val = np.asarray(val, dtype=float).ravel()

        ok = np.isfinite(pres) & np.isfinite(val)

        if ok.any():
            rows.append(pd.DataFrame({
                "LONG": float(lon),
                "LAT": float(lat),
                "PRES": pres[ok],
                "value": val[ok],
            }))

    if not rows:
        return pd.DataFrame(columns=["LONG", "LAT", "PRES", "value"])

    return pd.concat(rows, ignore_index=True)


def load_all_mat_points(matfiles, value_name):
    parts = [load_profiles_from_mat(f, value_name) for f in matfiles]
    return pd.concat(parts, ignore_index=True)


def make_statistics(df):
    edges = np.r_[0, 1, np.arange(3, np.ceil(df["PRES"].max()) + 2, 2)]

    df = df.copy()
    df["PRES_window"] = pd.cut(df["PRES"], bins=edges)

    stats = df.groupby("PRES_window", observed=False)["value"].describe()
    stats["midbin"] = stats.index.map(lambda x: x.mid).astype(float)

    return stats.reset_index()

def plot_selected_profiles(selected, area, out_path):
    """Plot accepted profile positions and the rectangular area box."""
    if selected.empty:
        return

    # one dot per profile/station, not one dot per depth level
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
    ax.set_title(f"{area['area']} accepted profiles, n={len(pts)}")
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

w_dir = Path(r"C:\Data\DMQC\UPDATE_test\\")

matfiles = [
    w_dir / "fmi_ctd_1601.mat",
    w_dir / "fmi_ctd_1602.mat",
    w_dir / "fmi_ctd_1501.mat",
    w_dir / "fmi_ctd_1502.mat",
]

areas = pd.read_csv(w_dir / "Climatology_boxes.txt")

variables = {
    "Practical_Salinity_dmnless": "SAL",
    "Temperature_degC": "TEMP",
}

for out_var_name, mat_var in variables.items():
    all_points = load_all_mat_points(matfiles, mat_var)

    for _, area in areas.iterrows():
        area_name = area["area"]

        selected = all_points[
            (all_points["LONG"] >= area["min_lon"]) &
            (all_points["LONG"] <= area["max_lon"]) &
            (all_points["LAT"] >= area["min_lat"]) &
            (all_points["LAT"] <= area["max_lat"])
        ]
        # sanity check to scrub obviously wrong data away
        if area_name == "BothSea" and out_var_name == "Practical_Salinity_dmnless":
            selected = selected[selected["value"] <= 20] #gludge to get rid of standard 35 psu value.
        if selected.empty:
            print(f"No data for {area_name} / {out_var_name}")
            continue

        stats = make_statistics(selected)

        outfile = w_dir / f"ICES_Statistics_{out_var_name}_by_depth_in_{area_name}.csv"
        stats.to_csv(outfile)  # keeps old unnamed index column
        mapfile = w_dir / f"Map_{out_var_name}_in_{area_name}.png"
        plot_selected_profiles(selected, area, mapfile)
        profilefile = w_dir / f"Profiles_{out_var_name}_in_{area_name}.png"
        plot_profile_cloud(selected, area_name, out_var_name, profilefile)

        print(outfile, selected.shape, stats.shape)