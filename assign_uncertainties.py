"""Ask for float-wide core uncertainty defaults and save YAML instructions.

Writes compact float-wide defaults, including future cycles.
No arguments uses config/local.yaml. This does not write D-files. Values use
native Argo units and apply to every profile index. Enter 'skip' to omit a
parameter; Ctrl-C cancels without saving. Existing measurements and QC are not
changed. Re-run dmqc_process.py write to apply saved instructions.
"""
import argparse
from pathlib import Path
import sys

from dmqc.settings import add_settings_argument, float_directory
from dmqc.instructions import CORE_PARAMETERS, parse_document, read_yaml, uncertainty_value
from dmqc.uncertainties import inspect_sources, write_defaults


def previous_defaults(path):
    if not path.exists():
        return {}, {}
    document = read_yaml(path)
    if document.get('checker') != 'default_uncertainties':
        raise ValueError(f'{path}: belongs to another checker; choose another output')
    entries = parse_document(document, path)
    values, reasons = {}, {}
    for p in CORE_PARAMETERS:
        matching = [s for s in entries if s.action == 'set_uncertainty' and p in s.target.parameters]
        scoped = [s for s in matching if s.scope == 'float']
        if scoped:
            matching = scoped
        unique = {s.value for s in matching}
        if len(unique) == 1:
            values[p] = unique.pop()
            reasons[p] = matching[0].reason
    return values, reasons


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_settings_argument(parser)
    parser.add_argument('directory', nargs='?', type=Path,
                        help='Float folder or R folder (default: local settings)')
    parser.add_argument('--output', type=Path, help='Default: instructions/uncertainties.yaml')
    parser.add_argument('--inspect', action='store_true', help='List existing estimates without prompting or saving')
    args = parser.parse_args(argv)
    try:
        args.directory = float_directory(args.directory, args.settings)
        r_dir, sources, stats = inspect_sources(args.directory)
        root = r_dir.parent if r_dir.name == 'R' else r_dir
        output = args.output or root / 'instructions' / 'uncertainties.yaml'
        if output.suffix.lower() not in ('.yaml', '.yml'):
            raise ValueError('Instruction output must be .yaml or .yml')
        print(f"{len(sources)} R-files, {sum(s['profiles'] for s in sources)} profile entries")
        for p, stat in stats.items():
            description = ('no usable existing estimates' if not stat['count'] else
                           f"{stat['count']} estimates, range {stat['minimum']:.7g} to {stat['maximum']:.7g}")
            print(f"{p} ({stat['units']}): {description}")
        if args.inspect:
            return 0
        saved, saved_reasons = previous_defaults(output)
        values, reasons = {}, {}
        print('These values replace existing uncertainties on eligible samples across all profiles.')
        print('No scientific values are assumed. Blank accepts a suggestion, or skips if none; skip removes a saved choice.')
        for p in stats:
            default = saved.get(p, stats[p]['default'])
            label = f' [{default:g}]' if default is not None else ''
            if default is not None:
                print(f"{p} suggestion from {'saved instructions' if p in saved else 'existing R-file estimates'}")
            while True:
                answer = input(f"{p} uncertainty ({stats[p]['units']}){label}: ").strip()
                if answer.lower() == 'skip' or (not answer and default is None):
                    break
                try:
                    value = uncertainty_value(default if not answer else float(answer))
                except ValueError as exc:
                    print(exc)
                    continue
                old_reason = saved_reasons.get(p, '') if value == saved.get(p) else ''
                reason = input(f"{p} source/justification{f' [{old_reason}]' if old_reason else ''}: ").strip()
                values[p] = value
                reasons[p] = reason or old_reason or 'Operator-selected default'
                break
        write_defaults(output, sources, values, reasons)
        print(f'Saved {values or "no defaults"} to {output}')
        return 0
    except (KeyboardInterrupt, EOFError):
        print('\nCancelled; instructions not saved.', file=sys.stderr)
        return 1
    except (ValueError, OSError, KeyError) as exc:
        print(f'Uncertainty setup failed: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
