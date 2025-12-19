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

def _cycle_from_filename(path):
    # expects R<id>_003.nc style
    name = Path(path).name
    parts = name.split("_")
    if len(parts) >= 2:
        return parts[1].split(".")[0]
    return "???"


def plot_cloud(ax, profiles, label, alpha=0.08, linewidth=0.7):
    for x, y in profiles:
        ax.plot(x, y, alpha=alpha, linewidth=linewidth, label=None)
    # add a dummy handle for legend
    ax.plot([], [], label=label)

def prepare_time_norm_and_labels(juld):
    """
    Given array of JULD values (float days since 1950-01-01, may contain nan),
    return: (juld_filled, norm, use_real_time, tick_positions, tick_labels)
    """
    juld = np.asarray(juld, dtype=float)

    if np.all(~np.isfinite(juld)):
        # fallback: file order proxy
        juld_filled = np.arange(len(juld), dtype=float)
        use_real_time = False
    else:
        use_real_time = True
        finite = juld[np.isfinite(juld)]
        jmin = float(np.min(finite))
        juld_filled = np.where(np.isfinite(juld), juld, jmin)

    norm = Normalize(vmin=float(np.min(juld_filled)), vmax=float(np.max(juld_filled)))

    tick_positions = None
    tick_labels = None
    if use_real_time:
        tick_positions = np.linspace(float(np.min(juld_filled)), float(np.max(juld_filled)), 6)
        tick_labels = [_juld_to_datetime(tt).strftime("%Y-%m-%d") for tt in tick_positions]

    return juld_filled, norm, use_real_time, tick_positions, tick_labels


def make_temp_psal_cloud_figure(
    rfiles,
    iprof=0,
    figsize=(12, 6),
    cmap_lines=None,
    alpha=0.1,
    linewidth=0.7,
):
    """
    Create (fig, (axT, axS), cbar) plotting TEMP and PSAL clouds from R-files only.
    Lines are colored by time (JULD) using the provided colormap.
    A single horizontal colorbar is placed under both panels.

    Returns:
        fig, (axT, axS), cbar
    """
    lines_temp = []
    lines_psal = []
    
    if cmap_lines is None:
        cmap_lines = LinearSegmentedColormap.from_list("psal_orange_green", ["orange", "green"])

    # --- Read times (JULD) for color mapping ---
    juld_raw = np.array([_read_juld(p, iprof=iprof) for p in rfiles], dtype=float)
    juld_filled, norm, use_real_time, tick_positions, tick_labels = prepare_time_norm_and_labels(juld_raw)

    # --- Build figure: TEMP left, PSAL right ---
    fig, (axT, axS) = plt.subplots(ncols=2, figsize=figsize, constrained_layout=True)

    n_T = 0
    n_S = 0

    for p, t in zip(rfiles, juld_filled):
        c = cmap_lines(norm(t))  # same colormap for both panels so one colorbar matches

        lt = None
        ls = None
        
        xT, yT = _read_profile_xy(p, "TEMP", y_var="PRES", iprof=iprof)
        if xT is not None:
            lt = axT.plot(xT, yT, color=c, linewidth=linewidth, alpha=alpha)[0]
            n_T += 1
        
        xS, yS = _read_profile_xy(p, "PSAL", y_var="PRES", iprof=iprof)
        if xS is not None:
            ls = axS.plot(xS, yS, color=c, linewidth=linewidth, alpha=alpha)[0]
            n_S += 1
        
        lines_temp.append(lt)
        lines_psal.append(ls)

    for ax in (axT, axS):
        ax.invert_yaxis()
        ax.grid(True, linewidth=0.4, alpha=0.5)

    axT.set_ylabel("Pressure (dbar)")
    axT.set_xlabel("Temperature (°C)")
    axT.set_title(f"TEMP cloud (R-files), iprof={iprof}  (n={n_T})")

    axS.set_xlabel("Salinity (PSU)")
    axS.set_title(f"PSAL cloud (R-files), iprof={iprof}  (n={n_S})")
    axS.set_ylabel("")

    # --- One shared colorbar (time), matching the plotted colours ---
    sm = ScalarMappable(norm=norm, cmap=cmap_lines)
    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=[axT, axS],
        orientation="horizontal",
        fraction=0.06,
        pad=0.08,
    )

    if use_real_time:
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)
        cbar.set_label("Profile time (JULD → date)")
    else:
        cbar.set_label("Profile order (JULD missing)")

    return fig, (axT, axS), cbar, lines_temp, lines_psal

def enable_profile_navigation(
    fig, axT, axS, rfiles, lines_temp, lines_psal,
    base_alpha=0.25, base_lw=0.7,
    sel_alpha=1.0, sel_lw=2.4
):
    """
    Adds up/down key navigation. Selected profile is drawn thicker & opaque.
    Keeps same color as original line.
    """

    n = len(rfiles)
    state = {"i": 0}

    def _apply_styles():
        i = state["i"]

        # reset all
        for lt, ls in zip(lines_temp, lines_psal):
            if lt is not None:
                lt.set_alpha(base_alpha)
                lt.set_linewidth(base_lw)
                lt.set_zorder(1)
            if ls is not None:
                ls.set_alpha(base_alpha)
                ls.set_linewidth(base_lw)
                ls.set_zorder(1)

        # highlight selected
        lt = lines_temp[i]
        ls = lines_psal[i]
        if lt is not None:
            lt.set_alpha(sel_alpha)
            lt.set_linewidth(sel_lw)
            lt.set_zorder(5)
        if ls is not None:
            ls.set_alpha(sel_alpha)
            ls.set_linewidth(sel_lw)
            ls.set_zorder(5)

        cyc = _cycle_from_filename(rfiles[i])
        axT.set_title(f"TEMP cloud (R-files) — selected #{i+1}/{n} (cycle {cyc})")
        axS.set_title(f"PSAL cloud (R-files) — selected #{i+1}/{n} (cycle {cyc})")

        fig.canvas.draw_idle()

    def _on_key(event):
        if event.key in ("down", "j"):
            state["i"] = (state["i"] - 1) % n
            _apply_styles()
        elif event.key in ("up", "k"):
            state["i"] = (state["i"] + 1) % n
            _apply_styles()

    fig.canvas.mpl_connect("key_press_event", _on_key)
    _apply_styles()
    return state  # returned so GUI can read current selection later



def main():
    ap = argparse.ArgumentParser(description="DMQC cloud plot: TEMP/PSAL vs PRES (R-files only)")
    ap.add_argument(
        "float_dir",
        nargs="?",
        type=str,
        default=r"C:\Data\ARGO_Dataa\DMQCprocessing\6903708\\",
        help="Path to float directory (contains R/)",
    )
    ap.add_argument("--iprof", type=int, default=0, help="Profile index inside each file (default: 0)")
    ap.add_argument("--save", type=str, default="", help="Output directory for PNGs (if empty: show interactively)")
    ap.add_argument("--dpi", type=int, default=180, help="PNG DPI if saving")
    args = ap.parse_args()

    float_dir = Path(args.float_dir)
    r_dir = float_dir / "R"
    rfiles = sorted(glob.glob(str(r_dir / "R*.nc")))

    if not rfiles:
        raise SystemExit(f"No R-files found in: {r_dir}")

    # colormap used both for lines and the shared colorbar
    cmap_psal = LinearSegmentedColormap.from_list("psal_orange_green", ["orange", "green"])

    fig, (axT, axS), cbar, lines_temp, lines_psal = make_temp_psal_cloud_figure(
        rfiles=rfiles,
        iprof=args.iprof,
        cmap_lines=cmap_psal,
    )
    
    enable_profile_navigation(fig, axT, axS, rfiles, lines_temp, lines_psal)

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
