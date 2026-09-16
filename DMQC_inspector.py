#!/usr/bin/env python3
import argparse
from pathlib import Path
from datetime import datetime, timedelta, timezone

from matplotlib.backend_bases import key_press_handler
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from dmqc.instructions import parse_document, read_yaml, sha256, write_profile_flags
ARGO_JULD_REF = datetime(1950, 1, 1)

def _get_profile_1d(var, iprof=0):
    """Return 1D array for iprof from a netCDF4 variable (supports 1D or 2D)."""
    if iprof < 0:
        return None
    if var.ndim == 2:
        if var.shape[0] <= iprof:
            return None
        return var[iprof, :].astype(np.float64)
    elif var.ndim == 1 and iprof == 0:
        return var[:].astype(np.float64)
    return None


def _mask_fill(arr, fill):
    if arr is None:
        return None
    a = np.ma.filled(np.ma.asarray(arr, dtype=np.float64), np.nan)
    valid = np.isfinite(a)
    if fill is not None:
        valid &= a != fill
    return np.where(valid, a, np.nan)


def _read_juld(nc_path, iprof=0):
    """Read JULD for profile iprof; returns float days since 1950-01-01 or np.nan."""
    with Dataset(nc_path, "r") as ds:
        if "JULD" not in ds.variables:
            return np.nan
        v = ds.variables["JULD"]
        # JULD has one timestamp per N_PROF entry, unlike depth-series data.
        if v.ndim == 0 or iprof < 0 or iprof >= v.shape[0]:
            return np.nan
        values = _mask_fill(v[iprof], getattr(v, "_FillValue", None)).ravel()
        if values.size != 1:
            return np.nan
        return float(values[0])


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


def _cycle_from_filename(path):
    # expects R<id>_003.nc style
    name = Path(path).name
    parts = name.split("_")
    if len(parts) >= 2:
        return parts[1].split(".")[0]
    return "???"


def available_profile_indices(rfiles):
    """Return indices present in at least one file, and each file's profile count."""
    counts = []
    for path in rfiles:
        with Dataset(path, "r") as ds:
            counts.append(len(ds.dimensions["N_PROF"]))
    return list(range(max(counts, default=0))), counts


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
    Plot TEMP and PSAL clouds from R-files only.
    Lines are colored by time (JULD) using the provided colormap.
    A single horizontal colorbar is placed under both panels.

    Returns:
        fig, (axT, axS), cbar, lines_temp, lines_psal
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

        xT, yT = _read_profile_xy(p, "TEMP", y_var="PRES", iprof=iprof)
        if xT is not None:
            n_T += 1
        else:
            xT, yT = [], []
        lt = axT.plot(xT, yT, color=c, linewidth=linewidth, alpha=alpha)[0]
        
        xS, yS = _read_profile_xy(p, "PSAL", y_var="PRES", iprof=iprof)
        if xS is not None:
            n_S += 1
        else:
            xS, yS = [], []
        ls = axS.plot(xS, yS, color=c, linewidth=linewidth, alpha=alpha)[0]
        
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
    fig, axT, axS, rfiles, lines_temp, lines_psal, cbar,
    base_alpha=0.25, base_lw=0.7,
    sel_alpha=1.0, sel_lw=2.4, iprof=0, instructions_path=None, source_hashes=None
):
    """
    Up/down selects a cycle; P cycles through available profile indices.
    Space toggles all indices in the selected file; F toggles only this index.
    Flagged profiles are red; Enter saves instructions without closing.
    Q saves and quits. Saved flags are restored when opening the inspector.
    """

    n = len(rfiles)
    indices, profile_counts = available_profile_indices(rfiles)
    if iprof not in indices:
        raise ValueError(f"Profile index {iprof} is unavailable; choose from {indices}")
    # Flags are keyed by source filename and profile index, not screen position.
    state = {"i": 0, "iprof": iprof, "flagged": set()}
    instructions_path = Path(instructions_path) if instructions_path is not None else (
        Path(rfiles[0]).parent.parent / "instructions" / "visual_inspector.yaml"
    )
    paths = {Path(path).name: Path(path) for path in rfiles}
    source_hashes = (dict(source_hashes) if source_hashes is not None
                     else {name: sha256(path) for name, path in paths.items()})
    counts_by_source = dict(zip(paths, profile_counts))
    for name, path in paths.items():
        if source_hashes.get(name) != sha256(path):
            raise ValueError(f"Source changed while loading: {name}; reopen the inspector")
    if instructions_path.exists():
        document = read_yaml(instructions_path)
        if document.get("checker") != "visual_inspector":
            raise ValueError(f"{instructions_path}: belongs to another checker")
        for instruction in parse_document(document, instructions_path):
            target = instruction.target
            if instruction.action != "flag":
                raise ValueError("Inspector report must contain only bad-profile flags")
            if target.source not in paths or target.profile_index >= counts_by_source[target.source]:
                raise ValueError(f"Saved target is unavailable: {target.source}, index {target.profile_index}")
            if instruction.source_sha256 is not None and instruction.source_sha256 != source_hashes[target.source]:
                raise ValueError(f"Saved source changed: {target.source}; review the saved report before continuing")
            state["flagged"].add((target.source, target.profile_index))
    saved_flags = set(state["flagged"])
    save_message = (f"Loaded {len(saved_flags)} saved flags" if instructions_path.exists()
                    else "No instructions saved yet")
    base_colors = [line.get_color() for line in lines_temp]
    dates = np.array([_read_juld(path, iprof=iprof) for path in rfiles])
    # Reuse already plotted data, and cache other indices after their first visit.
    cache = {iprof: (dates, [(lt.get_data(), ls.get_data())
                            for lt, ls in zip(lines_temp, lines_psal)])}
    date_title = fig.suptitle("")

    def _refresh_cloud():
        nonlocal dates, base_colors
        index = state["iprof"]
        if index not in cache:
            profile_data = []
            for path, count in zip(rfiles, profile_counts):
                if index < count:
                    if sha256(path) != source_hashes[Path(path).name]:
                        raise ValueError(f"Source changed: {Path(path).name}; reopen the inspector")
                    pairs = [_read_profile_xy(path, name, iprof=index)
                             for name in ("TEMP", "PSAL")]
                    profile_data.append([(x, y) if x is not None else ([], [])
                                         for x, y in pairs])
                else:
                    profile_data.append([([], []), ([], [])])
            cache[index] = (
                np.array([_read_juld(path, iprof=index) for path in rfiles]),
                profile_data,
            )
        dates, profile_data = cache[index]
        times, norm, use_dates, ticks, labels = prepare_time_norm_and_labels(dates)
        cbar.mappable.set_norm(norm)
        base_colors = [cbar.mappable.to_rgba(time) for time in times]
        for lt, ls, pairs, time in zip(lines_temp, lines_psal, profile_data, times):
            for line, (x, y) in zip((lt, ls), pairs):
                line.set_data(x, y)
                line.set_color(cbar.mappable.to_rgba(time))
        if use_dates:
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(labels)
            cbar.set_label("Profile time (JULD → date)")
        else:
            order_ticks = np.unique(np.linspace(times.min(), times.max(), min(n, 6)).astype(int))
            cbar.set_ticks(order_ticks)
            cbar.set_ticklabels([str(time) for time in order_ticks])
            cbar.set_label("Profile order (JULD missing)")
        for ax in (axT, axS):
            ax.relim()
            ax.autoscale(enable=True)
            ax.yaxis.set_inverted(True)

    def _apply_styles():
        i = state["i"]

        # Restore time colours when unflagged; red overrides colour only.
        for position, (lt, ls) in enumerate(zip(lines_temp, lines_psal)):
            key = (Path(rfiles[position]).name, state["iprof"])
            color = "#ff0000" if key in state["flagged"] else base_colors[position]
            if lt is not None:
                lt.set_color(color)
                lt.set_alpha(base_alpha)
                lt.set_linewidth(base_lw)
                lt.set_zorder(1)
            if ls is not None:
                ls.set_color(color)
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

        index = state["iprof"]
        date = (_juld_to_datetime(dates[i]).strftime("%Y-%m-%d %H:%M UTC")
                if np.isfinite(dates[i]) else "Date unavailable")
        availability = " — unavailable in this cycle" if index >= profile_counts[i] else ""
        key = (Path(rfiles[i]).name, index)
        status = "FLAGGED BAD" if key in state["flagged"] else "Unflagged"
        if availability:
            status = "Unavailable"
        date_title.set_text(
            f"Selected profile date: {date}\n"
            f"Profile index {index} of {indices}{availability} | {status} "
            f"| Total flagged: {len(state['flagged'])}\n"
            "Space: all indices | F: this index | P: change index | ↑/↓: cycle | Enter: save | Q: save & quit"
        )
        dirty = " — unsaved changes" if state["flagged"] != saved_flags else ""
        fig.supxlabel(f"{save_message}{dirty}", fontsize=9)
        cyc = _cycle_from_filename(rfiles[i])
        axT.set_title(f"TEMP cloud (R-files) — selected #{i+1}/{n} (cycle {cyc})")
        axS.set_title(f"PSAL cloud (R-files) — selected #{i+1}/{n} (cycle {cyc})")

        fig.canvas.draw_idle()

    def _save_flags():
        nonlocal saved_flags, save_message
        try:
            for name in {source for source, _ in state["flagged"]}:
                if sha256(paths[name]) != source_hashes[name]:
                    raise ValueError(f"Source changed: {name}; reopen the inspector")
            write_profile_flags(
                instructions_path, state["flagged"], source_hashes=source_hashes,
                reason="Whole profile rejected during visual inspection.",
                metadata={"created_utc": datetime.now(timezone.utc).isoformat()},
            )
        except (OSError, ValueError) as exc:
            save_message = f"Save failed: {exc}"
            print(save_message)
            _apply_styles()
            return False
        saved_flags = set(state["flagged"])
        save_message = f"Saved {len(saved_flags)} flags to {instructions_path}"
        print(save_message)
        _apply_styles()
        return True

    def _on_key(event):
        if event.key in ("enter", "return"):
            _save_flags()
        elif event.key in ("q", "Q"):
            if _save_flags():
                plt.close(fig)
        elif event.key in ("down", "j"):
            state["i"] = (state["i"] - 1) % n
            _apply_styles()
        elif event.key in ("up", "k"):
            state["i"] = (state["i"] + 1) % n
            _apply_styles()
        elif event.key in ("p", "P") and len(indices) > 1:
            position = indices.index(state["iprof"])
            state["iprof"] = indices[(position + 1) % len(indices)]
            _refresh_cloud()
            _apply_styles()
        elif event.key in (" ", "space", "f", "F"):
            count = profile_counts[state["i"]]
            if state["iprof"] < count:
                source = Path(rfiles[state["i"]]).name
                current = (source, state["iprof"])
                # Space makes every index match the displayed index's NEW state,
                # even when the indices previously had different flags.
                flag_bad = current not in state["flagged"]
                targets = range(count) if event.key in (" ", "space") else [state["iprof"]]
                for index in targets:
                    key = (source, index)
                    if flag_bad:
                        state["flagged"].add(key)
                    else:
                        state["flagged"].discard(key)
                _apply_styles()

    def _standard_keys(event):
        if event.key not in ("p", "P", " ", "space", "f", "F", "enter", "return", "q", "Q"):
            key_press_handler(event)

    manager = fig.canvas.manager
    if manager is not None:
        fig.canvas.mpl_disconnect(manager.key_press_handler_id)
        manager.key_press_handler_id = fig.canvas.mpl_connect("key_press_event", _standard_keys)
    fig.canvas.mpl_connect("key_press_event", _on_key)
    _apply_styles()
    return state  # returned so GUI can read current selection later



def main():
    ap = argparse.ArgumentParser(description="DMQC cloud plot: TEMP/PSAL vs PRES (R-files only)")
    ap.add_argument(
        "float_dir",
        nargs="?",
        type=str,
        default=r"/mnt/c/Data/ARGO_Dataa/DMQCprocessing/6903708/",
        help="Path to float directory (contains R/)",
    )
    ap.add_argument("--iprof", type=int, default=0, help="Profile index inside each file (default: 0)")
    ap.add_argument("--save", type=str, default="", help="Output directory for PNGs (if empty: show interactively)")
    ap.add_argument("--instructions", type=Path, help="Inspector YAML path (default: <float>/instructions/visual_inspector.yaml)")
    ap.add_argument("--dpi", type=int, default=180, help="PNG DPI if saving")
    args = ap.parse_args()

    float_dir = Path(args.float_dir)
    r_dir = float_dir / "R"
    rfiles = sorted(r_dir.glob("R*.nc"))

    if not rfiles:
        raise SystemExit(f"No R-files found in: {r_dir}")

    indices, _ = available_profile_indices(rfiles)
    if args.iprof not in indices:
        ap.error(f"Profile index {args.iprof} is unavailable; choose from {indices}")

    # Pin source versions before loading the plotted observations.
    source_hashes = {path.name: sha256(path) for path in rfiles}
    fig, (axT, axS), cbar, lines_temp, lines_psal = make_temp_psal_cloud_figure(
        rfiles=rfiles,
        iprof=args.iprof,
    )

    enable_profile_navigation(
        fig, axT, axS, rfiles, lines_temp, lines_psal, cbar, iprof=args.iprof,
        instructions_path=args.instructions or float_dir / "instructions" / "visual_inspector.yaml",
        source_hashes=source_hashes,
    )

    if args.save:
        outdir = Path(args.save)
        outdir.mkdir(parents=True, exist_ok=True)
        outpath = outdir / "cloud_TEMP_PSAL_R_only.png"
        fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {outpath}")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
