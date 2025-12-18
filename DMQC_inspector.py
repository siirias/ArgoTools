#!/usr/bin/env python3
import argparse
import glob
import os
from pathlib import Path
from datetime import datetime, timedelta

from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
ARGO_JULD_REF = datetime(1950, 1, 1)

def _get_profile_1d(var, iprof=0):
    """Return 1D array for iprof from a netCDF4 variable (supports 1D or 2D)."""
    if var.ndim == 2:
        if var.shape[0] <= iprof:
            return None
        return var[iprof, :].astype(np.float64)
    elif var.ndim == 1:
        return var[:].astype(np.float64)
    return None


def _mask_fill(arr, fill):
    if arr is None:
        return None
    a = np.array(arr, dtype=np.float64, copy=False)
    if fill is None:
        # Some Argo vars can be NaN-filled in practice, so keep this too:
        return np.where(np.isfinite(a), a, np.nan)
    return np.where(a != fill, a, np.nan)

def _read_juld(nc_path, iprof=0):
    """Read JULD for profile iprof; returns float days since 1950-01-01 or np.nan."""
    with Dataset(nc_path, "r") as ds:
        if "JULD" not in ds.variables:
            return np.nan
        v = ds.variables["JULD"]
        j = _get_profile_1d(v, iprof=iprof)
        if j is None or len(j) == 0:
            return np.nan
        # JULD is usually scalar per profile; but sometimes stored as length-1 arrays depending on templates
        val = float(np.ravel(j)[0])
        fill = getattr(v, "_FillValue", None)
        if fill is not None and val == float(fill):
            return np.nan
        if not np.isfinite(val):
            return np.nan
        return val
    
def _read_profile_xy(nc_path, x_var, y_var="PRES", iprof=0):
    """
    Read x and y as 1D arrays for profile iprof.
    Returns (x, y) with NaNs removed, or (None, None).
    """
    with Dataset(nc_path, "r") as ds:
        if x_var not in ds.variables or y_var not in ds.variables:
            return None, None

        vx = ds.variables[x_var]
        vy = ds.variables[y_var]

        x = _mask_fill(_get_profile_1d(vx, iprof=iprof), getattr(vx, "_FillValue", None))
        y = _mask_fill(_get_profile_1d(vy, iprof=iprof), getattr(vy, "_FillValue", None))

    if x is None or y is None:
        return None, None

    n = min(len(x), len(y))
    x = x[:n]
    y = y[:n]

    ok = np.isfinite(x) & np.isfinite(y)
    if not np.any(ok):
        return None, None
    return x[ok], y[ok]

def _juld_to_datetime(juld_days):
    return ARGO_JULD_REF + timedelta(days=float(juld_days))


def _collect_cloud(file_list, x_names, y_name="PRES", iprof=0):
    """Return list of (x, y) profiles found."""
    out = []
    used_names = {}
    for p in file_list:
        x, y, used = _read_profile_xy(p, x_names=x_names, y_name=y_name, iprof=iprof)
        if x is None:
            continue
        out.append((x, y))
        used_names[used] = used_names.get(used, 0) + 1
    return out, used_names


def plot_cloud(ax, profiles, label, alpha=0.08, linewidth=0.7):
    for x, y in profiles:
        ax.plot(x, y, alpha=alpha, linewidth=linewidth, label=None)
    # add a dummy handle for legend
    ax.plot([], [], label=label)


def main():
    ap = argparse.ArgumentParser(description="DMQC cloud plot: TEMP/PSAL vs PRES for all R and D profiles")
    ap.add_argument("float_dir", nargs='?', type=str, default = r"C:\Data\ARGO_Dataa\DMQCprocessing\6903708\\", help="Path to float directory (contains R/ and D/)")
    ap.add_argument("--iprof", type=int, default=0, help="Profile index inside each file (default: 0)")
    ap.add_argument("--save", type=str, default="", help="Output directory for PNGs (if empty: show interactively)")
    ap.add_argument("--dpi", type=int, default=180, help="PNG DPI if saving")
    args = ap.parse_args()

    float_dir = Path(args.float_dir)
    r_dir = float_dir / "R"
    d_dir = float_dir / "D"

    rfiles = sorted(glob.glob(str(r_dir / "R*.nc")))
    dfiles = sorted(glob.glob(str(d_dir / "D*.nc")))

    if not rfiles:
        raise SystemExit(f"No R-files found in: {r_dir}")
    if not dfiles:
        print(f"Warning: no D-files found in: {d_dir} (will plot only R)")

    # --- Read times (JULD) for color mapping ---
    juld = np.array([_read_juld(p, iprof=args.iprof) for p in rfiles], dtype=float)

    # Fallback if JULD missing everywhere: use file order as a proxy "time"
    if np.all(~np.isfinite(juld)):
        juld = np.arange(len(rfiles), dtype=float)
        use_real_time = False
    else:
        use_real_time = True
        # If some files have missing JULD, shove them to the beginning (oldest)
        finite = juld[np.isfinite(juld)]
        jmin = float(np.min(finite))
        juld = np.where(np.isfinite(juld), juld, jmin)

    norm = Normalize(vmin=float(np.min(juld)), vmax=float(np.max(juld)))

    # Custom colormaps
    cmap_temp = LinearSegmentedColormap.from_list("temp_black_red", ["black", "red"])
    cmap_psal = LinearSegmentedColormap.from_list("psal_orange_green", ["orange", "green"])

    # --- Build figure: TEMP left, PSAL right ---
    fig, (axT, axS) = plt.subplots(ncols=2, figsize=(12, 6), constrained_layout=True)

    n_T = 0
    n_S = 0

    for p, t in zip(rfiles, juld):
        cS = cmap_psal(norm(t))
        cT = cS #Identical colormap, for comparison
        xT, yT = _read_profile_xy(p, "TEMP", y_var="PRES", iprof=args.iprof)
        if xT is not None:
            axT.plot(xT, yT, color=cT, linewidth=0.7, alpha=0.25)
            n_T += 1

        xS, yS = _read_profile_xy(p, "PSAL", y_var="PRES", iprof=args.iprof)
        if xS is not None:
            axS.plot(xS, yS, color=cS, linewidth=0.7, alpha=0.25)
            n_S += 1

    for ax in (axT, axS):
        ax.invert_yaxis()
        ax.set_ylabel("Pressure (dbar)")
        ax.grid(True, linewidth=0.4, alpha=0.5)

    axT.set_xlabel("Temperature (°C)")
    axT.set_title(f"TEMP cloud (R-files), iprof={args.iprof}  (n={n_T})")

    axS.set_xlabel("Salinity (PSU)")
    axS.set_title(f"PSAL cloud (R-files), iprof={args.iprof}  (n={n_S})")
    axS.set_ylabel("")  # keep only left ylabel for cleanliness

    # --- One shared colorbar (time) ---
    # Use a neutral mappable for the colorbar scale; the plotted lines use different colormaps.
    sm = ScalarMappable(norm=norm, cmap=plt.cm.viridis)
    sm.set_array([])

    # --- One shared colorbar (time), matching the plotted colours ---
    sm = ScalarMappable(norm=norm, cmap=cmap_psal)
    sm.set_array([])
    
    cbar = fig.colorbar(
        sm,
        ax=[axT, axS],
        orientation="horizontal",
        fraction=0.06,
        pad=0.08,
    )
    
    if use_real_time:
        ticks = np.linspace(float(np.min(juld)), float(np.max(juld)), 6)
        cbar.set_ticks(ticks)
        labels = [_juld_to_datetime(tt).strftime("%Y-%m-%d") for tt in ticks]
        cbar.set_ticklabels(labels)
        cbar.set_label("Profile time (JULD → date)")
    else:
        cbar.set_label("Profile order (JULD missing)")
    if args.save:
        outdir = Path(args.save)
        outdir.mkdir(parents=True, exist_ok=True)
        outpath = outdir / "cloud_TEMP_PSAL_R_only.png"
        fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {outpath}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
