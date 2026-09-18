#!/usr/bin/env python3
import argparse
import textwrap
from pathlib import Path
from datetime import datetime, timedelta, timezone

from matplotlib.backend_bases import key_press_handler
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from dmqc.instructions import (parse_document, read_yaml, sha256, write_instructions,
                               load_instructions, Target, Instruction, Decision)
from dmqc.download import file_identity
from dmqc.profiles import profile_parameters
from dmqc.settings import add_settings_argument, float_directory
ARGO_JULD_REF = datetime(1950, 1, 1)

def _mask_fill(arr, fill):
    if arr is None:
        return None
    a = np.ma.filled(np.ma.asarray(arr, dtype=np.float64), np.nan)
    valid = np.isfinite(a)
    if fill is not None:
        valid &= a != fill
    return np.where(valid, a, np.nan)


def _juld_to_datetime(juld_days):
    return ARGO_JULD_REF + timedelta(days=float(juld_days))


def _cycle_from_filename(path):
    # expects R<id>_003.nc style
    name = Path(path).name
    parts = name.split("_")
    if len(parts) >= 2:
        return parts[1].split(".")[0]
    return "???"


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
        tick_positions = np.unique(np.linspace(float(np.min(juld_filled)), float(np.max(juld_filled)), 6))
        tick_labels = [_juld_to_datetime(tt).strftime("%Y-%m-%d") for tt in tick_positions]

    return juld_filled, norm, use_real_time, tick_positions, tick_labels


class ProfileCloud:
    """Metadata-driven raw measurement panels, cached by profile index."""
    AUXILIARY = {'PRES', 'MTIME', 'NB_SAMPLE_CTD'}

    def __init__(self, rfiles, iprof=0, source_hashes=None):
        self.rfiles = [Path(p) for p in rfiles]
        self.hashes = dict(source_hashes) if source_hashes is not None else {
            p.name: sha256(p) for p in self.rfiles}
        self.inventory, self.labels, self.counts, self.cache = [], {}, [], {}
        self.notices = set()
        pressure_units = set()
        for path in self.rfiles:
            with Dataset(path) as ds:
                count = len(ds.dimensions['N_PROF'])
                self.counts.append(count)
                profiles = []
                for ip in range(count):
                    names = profile_parameters(ds, ip)
                    measurements = {}
                    if 'PRES' not in names or 'PRES' not in ds.variables:
                        self.notices.add(f'{path.name}: profile {ip} has no pressure coordinate')
                        profiles.append(measurements)
                        continue
                    pressure = ds['PRES']
                    if pressure.dimensions != ('N_PROF', 'N_LEVELS'):
                        raise ValueError(f'{path.name}: unsupported PRES dimensions {pressure.dimensions}')
                    unit = str(getattr(pressure, 'units', '')).strip()
                    pressure_units.add('dbar' if unit in ('decibar', 'decibars', 'dbar') else unit)
                    for name in names:
                        if name in self.AUXILIARY or name.endswith(('_QC', '_ERROR', '_ADJUSTED', '_MED', '_STD')):
                            continue
                        if name not in ds.variables:
                            self.notices.add(f'{path.name}: declared {name} has no variable; not plotted')
                            continue
                        var = ds[name]
                        if var.dimensions != pressure.dimensions or var.shape != pressure.shape or var.dtype.kind not in 'fiu':
                            self.notices.add(f'{path.name}: {name} does not share the pressure sampling dimensions; not plotted')
                            continue
                        units = str(getattr(var, 'units', '')).strip()
                        key = (name, units)
                        measurements[key] = name
                        self.labels.setdefault(key, str(getattr(var, 'long_name', name)).strip() or name)
                    profiles.append(measurements)
                self.inventory.append(profiles)
            if sha256(path) != self.hashes[path.name]:
                raise ValueError(f'Source changed while loading: {path.name}')
        if len(pressure_units) > 1:
            raise ValueError(f'Inconsistent pressure units across files: {sorted(pressure_units)}; cannot share a pressure axis')
        self.pressure_unit = next(iter(pressure_units), '') or 'units unavailable'
        self.indices = list(range(max(self.counts, default=0)))
        if iprof not in self.indices:
            raise ValueError(f'Profile index {iprof} is unavailable; choose from {self.indices}')
        self.fig = plt.figure(figsize=(12, 6))
        self.axes, self.lines, self.highlights, self.keys = [], [], [], None
        self.cbar = None
        self.cmap = LinearSegmentedColormap.from_list('profile_time', ['orange', 'green'])
        self.show_index(iprof)
        for notice in sorted(self.notices):
            print(f'Viewer: {notice}')

    def _load_index(self, index):
        keys = sorted({key for profiles in self.inventory if index < len(profiles)
                       for key in profiles[index]},
                      key=lambda key: (0 if key[0] == 'TEMP' else 1 if key[0] == 'PSAL' else 2, *key))
        dates, data = [], []
        for path, count in zip(self.rfiles, self.counts):
            pairs = {}
            date = np.nan
            if index < count:
                if sha256(path) != self.hashes[path.name]:
                    raise ValueError(f'Source changed: {path.name}; reopen the inspector')
                with Dataset(path) as ds:
                    if 'JULD' in ds.variables:
                        v = ds['JULD']
                        if v.dimensions == ('N_PROF',):
                            date = float(_mask_fill(v[index], getattr(v, '_FillValue', None)))
                    declared = profile_parameters(ds, index)
                    if 'PRES' in declared and 'PRES' in ds.variables:
                        v = ds['PRES']
                        pressure = _mask_fill(v[index], getattr(v, '_FillValue', None))
                        for key in keys:
                            name, units = key
                            if name not in declared or name not in ds.variables:
                                continue
                            v = ds[name]
                            if (v.dimensions != ('N_PROF', 'N_LEVELS') or v.shape != ds['PRES'].shape
                                    or v.dtype.kind not in 'fiu' or str(getattr(v, 'units', '')).strip() != units):
                                continue
                            values = _mask_fill(v[index], getattr(v, '_FillValue', None))
                            valid = np.isfinite(values) & np.isfinite(pressure)
                            pairs[key] = (values[valid], pressure[valid])
                if sha256(path) != self.hashes[path.name]:
                    raise ValueError(f'Source changed: {path.name}; reopen the inspector')
            dates.append(date)
            data.append(pairs)
        self.cache[index] = keys, np.asarray(dates), data

    def show_index(self, index):
        if index not in self.cache:
            self._load_index(index)
        keys, self.dates, self.data = self.cache[index]
        # An empty panel explains pressure-only profiles without inventing a curve.
        keys = keys or [(None, '')]
        rebuilt = keys != self.keys
        if rebuilt:
            for ax in self.axes:
                ax.remove()
            self.keys = keys
            columns = min(3, len(keys))
            rows = (len(keys) + columns - 1) // columns
            height = 6 + 3 * (rows - 1)
            self.fig.set_size_inches(12 if columns < 3 else 15, height, forward=True)
            grid = self.fig.add_gridspec(rows, columns, left=.07, right=.98,
                                         bottom=2.45 / height, top=1 - 1.1 / height,
                                         wspace=.30, hspace=.60)
            self.axes, self.lines, self.highlights = [], [], []
            for j, key in enumerate(keys):
                ax = self.fig.add_subplot(grid[j // columns, j % columns],
                                          sharey=self.axes[0] if self.axes else None)
                self.axes.append(ax)
                self.lines.append([ax.plot([], [], linewidth=.7, alpha=.25, zorder=1)[0] for _ in self.rfiles])
                self.highlights.append(ax.plot([], [], linewidth=2.4, alpha=1, zorder=5)[0])
                ax.grid(True, linewidth=.4, alpha=.5)
                ax.set_ylabel(f'Pressure ({self.pressure_unit})')
                if key[0] is not None:
                    label = textwrap.fill(self.labels[key], width=38)
                    ax.set_xlabel(f'{label}\n({key[1] or "units unavailable"})', fontsize=9)
                else:
                    ax.set_xlabel('Pressure only or no compatible measurement arrays')
            if self.cbar is not None:
                self.cbar.ax.set_position([.22, 1.35 / height, .56, .15 / height])
        times, norm, use_dates, ticks, labels = prepare_time_norm_and_labels(self.dates)
        self.colors = [self.cmap(norm(t)) for t in times]
        if self.cbar is None:
            sm = ScalarMappable(norm=norm, cmap=self.cmap)
            self.cbar = self.fig.colorbar(sm, cax=self.fig.add_axes([.22, 1.35 / self.fig.get_figheight(), .56,
                                                                  .15 / self.fig.get_figheight()]), orientation='horizontal')
        self.cbar.mappable.set_norm(norm)
        if use_dates:
            self.cbar.set_ticks(ticks)
            self.cbar.set_ticklabels(labels)
            self.cbar.set_label('Profile time (JULD → date)')
        else:
            ticks = np.unique(np.linspace(times.min(), times.max(), min(len(times), 6)).astype(int))
            self.cbar.set_ticks(ticks)
            self.cbar.set_ticklabels([str(t) for t in ticks])
            self.cbar.set_label('Profile order (JULD missing)')
        pressures = []
        for ax, key, lines, overlay in zip(self.axes, self.keys, self.lines, self.highlights):
            overlay.set_data([], [])
            for line, pairs, color in zip(lines, self.data, self.colors):
                x, y = pairs.get(key, ([], []))
                line.set_data(x, y)
                line.set_color(color)
                if len(y):
                    pressures.append(y)
            ax.relim()
            ax.autoscale(enable=True, axis='x')
        if pressures:
            low = min(float(y.min()) for y in pressures)
            high = max(float(y.max()) for y in pressures)
            margin = max((high - low) * .05, .1)
            self.axes[0].set_ylim(high + margin, low - margin)
        else:
            self.axes[0].set_ylim(1, 0)
        return rebuilt


class _NavigationRenderer:
    """Cache the static cloud; repaint highlights and text on navigation.

    Artists stay non-animated so normal savefig/toolbar exports include them.
    Only rebuilding the background needs a second draw without the overlays.
    """
    def __init__(self, fig, artists, axes):
        self.fig = fig
        self.artists = artists
        self.background = None
        self.capturing = False
        self.supports_blit = fig.canvas.supports_blit
        fig.canvas.mpl_connect('draw_event', self._on_draw)
        fig.canvas.mpl_connect('resize_event', self.invalidate)
        self.watch_axes(axes)

    def watch_axes(self, axes):
        for ax in axes:
            ax.callbacks.connect('xlim_changed', self.invalidate)
            ax.callbacks.connect('ylim_changed', self.invalidate)

    def invalidate(self, *_):
        self.background = None

    def _on_draw(self, event):
        if self.capturing or not self.supports_blit:
            return
        self.invalidate()
        if event.canvas.is_saving():
            return
        visible = [artist.get_visible() for artist in self.artists]
        self.capturing = True
        try:
            for artist in self.artists:
                artist.set_visible(False)
            event.canvas.draw()
            self.background = event.canvas.copy_from_bbox(self.fig.bbox)
        finally:
            for artist, was_visible in zip(self.artists, visible):
                artist.set_visible(was_visible)
            self.capturing = False
        self.paint()

    def paint(self, full=False):
        if full:
            self.invalidate()
        if not self.supports_blit or self.background is None:
            self.fig.canvas.draw_idle()
            return
        self.fig.canvas.restore_region(self.background)
        for artist in self.artists:
            if artist.get_visible():
                self.fig.draw_artist(artist)
        self.fig.canvas.blit(self.fig.bbox)


def enable_profile_navigation(
    cloud,
    base_alpha=0.25, base_lw=0.7,
    sel_alpha=1.0, sel_lw=2.4, iprof=0, instructions_path=None, source_hashes=None,
    instructions_dir=None
):
    """
    Up/down selects a cycle; P cycles through available profile indices.
    Space toggles all indices in the selected file; F toggles only this index.
    Other checkers are loaded first; visual accept/flag decisions have priority 100.
    R clears this index; Shift+R clears all indices in the cycle.
    Rejected profiles are red; Enter saves, Q saves and quits.
    """

    fig, rfiles = cloud.fig, cloud.rfiles
    n = len(rfiles)
    indices, profile_counts = cloud.indices, cloud.counts
    if iprof not in indices:
        raise ValueError(f"Profile index {iprof} is unavailable; choose from {indices}")
    # Flags are keyed by source filename and profile index, not screen position.
    state = {"i": 0, "iprof": iprof, "flagged": set(), "overrides": {}}
    instructions_path = Path(instructions_path) if instructions_path is not None else (
        Path(rfiles[0]).parent.parent / "instructions" / "visual_inspector.yaml"
    )
    paths = {Path(path).name: Path(path) for path in rfiles}
    source_hashes = (dict(source_hashes) if source_hashes is not None
                     else cloud.hashes)
    counts_by_source = dict(zip(paths, profile_counts))
    for name, path in paths.items():
        if source_hashes.get(name) != sha256(path):
            raise ValueError(f"Source changed while loading: {name}; reopen the inspector")
    instructions_dir = Path(instructions_dir) if instructions_dir is not None else (
        Path(rfiles[0]).parent.parent / 'instructions' if Path(rfiles[0]).parent.name == 'R'
        else instructions_path.parent)

    def _other_report_hashes():
        return {path.resolve(): sha256(path) for pattern in ('*.yaml', '*.yml')
                for path in instructions_dir.rglob(pattern)
                if path.resolve() != instructions_path.resolve()}

    other_hashes = _other_report_hashes()
    others = load_instructions(instructions_dir, Path(rfiles[0]).parent,
                              file_identity(Path(rfiles[0]).name)[0], exclude_paths=[instructions_path])
    external = {}
    for instruction in others:
        if instruction.action not in ('flag', 'accept', 'no_finding'):
            continue
        target = instruction.target
        if target.source not in paths or target.profile_index >= counts_by_source[target.source]:
            raise ValueError(f'Instruction target unavailable: {target}')
        if instruction.source_sha256 is not None and instruction.source_sha256 != source_hashes[target.source]:
            raise ValueError(f'Source changed since {instruction.checker}: {target.source}')
        external.setdefault((target.source, target.profile_index), []).append(instruction)
    if _other_report_hashes() != other_hashes:
        raise ValueError('Instructions changed while loading; reopen the inspector')
    if instructions_path.exists():
        document = read_yaml(instructions_path)
        if document.get("checker") != "visual_inspector":
            raise ValueError(f"{instructions_path}: belongs to another checker")
        for instruction in parse_document(document, instructions_path):
            target = instruction.target
            if instruction.action not in ("flag", "accept"):
                raise ValueError("Inspector report must contain only accept/flag decisions")
            if target.source not in paths or target.profile_index >= counts_by_source[target.source]:
                raise ValueError(f"Saved target is unavailable: {target.source}, index {target.profile_index}")
            if instruction.source_sha256 is not None and instruction.source_sha256 != source_hashes[target.source]:
                raise ValueError(f"Saved source changed: {target.source}; review the saved report before continuing")
            key = (target.source, target.profile_index)
            if key in state['overrides']:
                raise ValueError(f'Duplicate visual decision: {key}')
            state['overrides'][key] = instruction
    saved_overrides = dict(state['overrides'])
    save_message = (f"Loaded {len(saved_overrides)} saved decisions" if instructions_path.exists()
                    else "No visual decisions saved yet")

    def _decision(key):
        suggestions = list(external.get(key, ()))
        if key in state['overrides']:
            suggestions.append(state['overrides'][key])
        return Decision(Target(*key), source_hashes[key[0]], tuple(suggestions))

    date_title = fig.suptitle("", y=0.99, fontsize=11)
    status_text = fig.supxlabel('', y=0.015, fontsize=9)
    renderer = _NavigationRenderer(fig, [], cloud.axes)
    fig._dmqc_renderer = renderer
    fig._dmqc_highlights = cloud.highlights
    fig._dmqc_cloud = cloud

    def _dynamic_artists():
        renderer.artists = [*cloud.highlights, date_title, status_text, *(ax.title for ax in cloud.axes)]
        fig._dmqc_highlights = cloud.highlights

    _dynamic_artists()
    source_names = [Path(path).name for path in rfiles]
    resolved = {}

    def _resolve_flags():
        resolved.clear()
        for key in set(external) | set(state['overrides']):
            resolved[key] = _decision(key)
        state['flagged'] = {key for key, decision in resolved.items() if decision.rejected}

    def _selected_decision(key):
        return resolved.get(key) or Decision(Target(*key), source_hashes[key[0]], ())


    def _refresh_cloud():
        renderer.invalidate()
        # Resizing during a layout rebuild may emit draw events on GUI backends.
        renderer.capturing = True
        try:
            rebuilt = cloud.show_index(state['iprof'])
            if rebuilt:
                renderer.watch_axes(cloud.axes)
            _dynamic_artists()
        finally:
            renderer.capturing = False

    def _apply_styles(full=False, decisions_changed=False):
        if decisions_changed:
            _resolve_flags()
        i = state['i']
        for key, lines, overlay in zip(cloud.keys, cloud.lines, cloud.highlights):
            if full:
                for position, line in enumerate(lines):
                    target = (source_names[position], state['iprof'])
                    # Auxiliary sensor panels are display-only, not marked bad by CTD decisions.
                    color = '#ff0000' if key[0] in ('TEMP', 'PSAL') and target in state['flagged'] else cloud.colors[position]
                    line.set_color(color)
                    line.set_alpha(base_alpha)
                    line.set_linewidth(base_lw)
                    line.set_zorder(1)
            overlay.set_data(*lines[i].get_data())
            overlay.set_color(lines[i].get_color())
            overlay.set_alpha(sel_alpha)
            overlay.set_linewidth(sel_lw)

        index = state["iprof"]
        date = (_juld_to_datetime(cloud.dates[i]).strftime("%Y-%m-%d %H:%M UTC")
                if np.isfinite(cloud.dates[i]) else "Date unavailable")
        availability = " — unavailable in this cycle" if index >= profile_counts[i] else ""
        key = (Path(rfiles[i]).name, index)
        status = "FLAGGED BAD" if key in state["flagged"] else "Unflagged"
        winners = _selected_decision(key).winning_suggestions if not availability else ()
        if winners:
            status += f' | Effective priority {winners[0].priority}'
        override = state['overrides'].get(key)
        if override is not None:
            status += f" | Visual: {override.action} (priority {override.priority})"
        elif external.get(key):
            status += " | From other checkers"
        if availability:
            status = "Unavailable"
        date_title.set_text(
            f"Selected profile date: {date}\n"
            f"Profile index {index} of {indices}{availability} | {status} "
            f"| Total flagged: {len(state['flagged'])}\n"
            "Space/F: core QC all/this | r/R: reset this/all | P: index | ↑/↓: cycle | Enter: save | Q: save & quit"
        )
        dirty = " — unsaved changes" if state["overrides"] != saved_overrides else ""
        details = [f'{s.checker}: {s.action}, priority {s.priority}: {s.reason or "No reason supplied"}'
                   for s in _selected_decision(key).suggestions] if not availability else []
        description = '\n'.join(textwrap.fill(line, width=130) for line in details)
        footer = f"{save_message}{dirty}\n{description}"
        status_text.set_text(footer)
        # Fit all reasons inside the reserved footer band without moving the axes.
        lines = max(1, len(footer.splitlines()))
        status_text.set_y(.08 / fig.get_figheight())
        status_text.set_fontsize(min(9, 72 * .72 / (lines * 1.2)))
        cyc = _cycle_from_filename(rfiles[i])
        for ax, panel in zip(cloud.axes, cloud.keys):
            name = panel[0]
            note = ''
            if name is None:
                title = 'No plottable measurement parameters'
            else:
                title = name if name in ('TEMP', 'PSAL') else f'{name} (view only)'
                if panel not in cloud.data[i]:
                    note = '\nNot measured in this profile'
                elif len(cloud.data[i][panel][0]) == 0:
                    note = '\nNo usable samples'
            ax.set_title(f'{title} — cycle {cyc}{note}', y=1.0, fontsize=10)

        renderer.paint(full=full)

    def _save_flags():
        nonlocal saved_overrides, save_message
        try:
            if _other_report_hashes() != other_hashes:
                raise ValueError('Other instructions changed; reopen the inspector before saving')
            for name in {source for source, _ in state["overrides"]}:
                if sha256(paths[name]) != source_hashes[name]:
                    raise ValueError(f"Source changed: {name}; reopen the inspector")
            entries = []
            for key, instruction in sorted(state['overrides'].items()):
                entry = {'target': {'source': key[0], 'profile_index': key[1], 'selection': 'whole_profile'},
                         'action': instruction.action, 'priority': instruction.priority,
                         'reason': instruction.reason, 'source_sha256': source_hashes[key[0]]}
                if instruction.action == 'flag':
                    entry['flag'] = '4'
                entries.append(entry)
            write_instructions(instructions_path, 'visual_inspector', entries,
                               metadata={"created_utc": datetime.now(timezone.utc).isoformat()})
        except (OSError, ValueError) as exc:
            save_message = f"Save failed: {exc}"
            print(save_message)
            _apply_styles()
            return False
        saved_overrides = dict(state["overrides"])
        save_message = f"Saved {len(saved_overrides)} decisions to {instructions_path}"
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
            _apply_styles(full=True)
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
                    action = 'flag' if flag_bad else 'accept'
                    state['overrides'][key] = Instruction(
                        Target(source, index), action, 'visual_inspector',
                        'Whole profile rejected during visual inspection.' if flag_bad else
                        'Whole-profile rejection overridden during visual inspection; source QC retained.',
                        str(instructions_path), flag='4' if flag_bad else None,
                        source_sha256=source_hashes[source], priority=100)
                _apply_styles(full=True, decisions_changed=True)

        elif event.key in ('r', 'R'):
            source = Path(rfiles[state['i']]).name
            targets = range(profile_counts[state['i']]) if event.key == 'R' else [state['iprof']]
            for index in targets:
                state['overrides'].pop((source, index), None)
            _apply_styles(full=True, decisions_changed=True)

    def _standard_keys(event):
        if event.key not in ("p", "P", " ", "space", "f", "F", "enter", "return", "q", "Q", "r", "R"):
            key_press_handler(event)

    manager = fig.canvas.manager
    if manager is not None:
        fig.canvas.mpl_disconnect(manager.key_press_handler_id)
        manager.key_press_handler_id = fig.canvas.mpl_connect("key_press_event", _standard_keys)
    fig.canvas.mpl_connect("key_press_event", _on_key)
    _apply_styles(full=True, decisions_changed=True)
    return state  # returned so GUI can read current selection later



def main(argv=None):
    ap = argparse.ArgumentParser(description="DMQC profile clouds: available measurements vs pressure (R-files only)")
    add_settings_argument(ap)
    ap.add_argument(
        "float_dir",
        nargs="?",
        type=Path,
        help="Path to float directory containing R/ (default: local settings)",
    )
    ap.add_argument("--iprof", type=int, default=0, help="Profile index inside each file (default: 0)")
    ap.add_argument("--save", type=str, default="", help="Output directory for PNGs (if empty: show interactively)")
    ap.add_argument("--instructions", type=Path, help="Inspector YAML path (default: <float>/instructions/visual_inspector.yaml)")
    ap.add_argument("--instructions-dir", type=Path, help="Directory containing other checker instructions")
    ap.add_argument("--dpi", type=int, default=180, help="PNG DPI if saving")
    args = ap.parse_args(argv)

    try:
        float_dir = float_directory(args.float_dir, args.settings)
    except (ValueError, OSError) as exc:
        ap.error(str(exc))
    r_dir = float_dir / "R"
    rfiles = sorted(r_dir.glob("R*.nc"))

    if not rfiles:
        raise SystemExit(f"No R-files found in: {r_dir}")

    # Pin source versions before loading the plotted observations.
    source_hashes = {path.name: sha256(path) for path in rfiles}
    try:
        cloud = ProfileCloud(rfiles, iprof=args.iprof, source_hashes=source_hashes)
        enable_profile_navigation(
            cloud, iprof=args.iprof,
            instructions_path=args.instructions or float_dir / "instructions" / "visual_inspector.yaml",
            source_hashes=source_hashes, instructions_dir=args.instructions_dir,
        )
    except (ValueError, OSError) as exc:
        ap.error(str(exc))
    fig = cloud.fig

    if args.save:
        outdir = Path(args.save)
        outdir.mkdir(parents=True, exist_ok=True)
        outpath = outdir / "cloud_profiles_R_only.png"
        fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {outpath}")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
