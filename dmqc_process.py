"""Argo DMQC workflow: download sources, combine checker reports, write D-files.

Run without arguments for float 6903708 in the historical work directory.
Use --help for individual stages. Scientific checkers run separately.
"""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import sys

from dmqc.download import download_r_files
from dmqc.combine import combine_instructions
from dmqc.instructions import load_instructions, read_yaml
from dmqc.writer import DEFAULT_META, write_d_files


DEFAULT_FLOAT = '6903708'
DEFAULT_WORK_DIR = Path(r'C:\Data\ARGO_Dataa\DMQCprocessing' if os.name == 'nt'
                        else '/mnt/c/Data/ARGO_Dataa/DMQCprocessing')


def parse_cycles(value):
    """Accept comma-separated cycle numbers and inclusive ranges, e.g. 1,3-5."""
    cycles = set()
    try:
        for part in value.split(','):
            bounds = part.split('-')
            if len(bounds) == 1:
                first = last = int(bounds[0])
            elif len(bounds) == 2:
                first, last = map(int, bounds)
            else:
                raise ValueError
            if first < 0 or last < first:
                raise ValueError
            cycles.update(range(first, last + 1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Use cycles such as 1,3-5') from exc
    return cycles


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('stage', nargs='?', default='all', choices=('download', 'combine', 'write', 'all'),
                        help='Default: download then combine/write available decisions')
    result.add_argument('--float', dest='float_id', default=DEFAULT_FLOAT)
    result.add_argument('--work-dir', type=Path, default=DEFAULT_WORK_DIR,
                        help='Parent of the float directory (default: %(default)s)')
    result.add_argument('--dac', default='coriolis')
    result.add_argument('--cycles', type=parse_cycles, help='Subset, e.g. 1,3-5')
    result.add_argument('--instructions-dir', type=Path,
                        help='Partner YAML directory; default: <float>/instructions')
    result.add_argument('--dry-run', action='store_true', help='Read/validate and print the plan; write nothing')
    result.add_argument('--overwrite', action='store_true', help='Replace existing D-files and their reports')
    return result


def run(args):
    if not args.float_id.isdigit():
        raise ValueError('--float must contain only digits')
    float_dir = args.work_dir / args.float_id
    r_dir, d_dir = float_dir / 'R', float_dir / 'D'
    log = sys.stderr if args.stage == 'combine' else sys.stdout
    print(f'Float {args.float_id}; stage={args.stage}; directory={float_dir}', file=log)
    if args.stage in ('download', 'all'):
        download_r_files(args.float_id, r_dir, dac=args.dac, cycles=args.cycles, dry_run=args.dry_run)
    if args.stage == 'download':
        return []
    instruction_dir = args.instructions_dir or float_dir / 'instructions'
    if args.instructions_dir is not None and not instruction_dir.is_dir():
        raise ValueError(f'Instructions directory does not exist: {instruction_dir}')
    instructions = load_instructions(instruction_dir, r_dir, args.float_id, args.cycles,
                                     legacy_dir=float_dir / 'cycles')
    decisions = combine_instructions(instructions, r_dir)
    summary = {'float_id': args.float_id, 'instruction_count': len(instructions),
               'bad_profile_count': len(decisions), 'decisions': [d.as_dict() for d in decisions],
               'no_findings': [asdict(i) for i in instructions if i.action == 'no_finding']}
    print(f'{len(instructions)} checker reports; {len(decisions)} rejected profiles', file=log)
    if args.stage == 'combine':
        # A reviewable plan on stdout, without mixing generated plans into checker inputs.
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return summary
    if not decisions:
        print(f'No bad-profile instructions; no D-files written. Add checker YAMLs in {instruction_dir}')
        return []
    meta = dict(DEFAULT_META)
    meta_path = float_dir / 'meta.yaml'
    if meta_path.exists():
        # Old processing defaults must not change decisions, errors, or mode fields.
        stored = read_yaml(meta_path)
        meta.update({k: stored[k] for k in DEFAULT_META if k in stored})
    for key, value in meta.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'meta.yaml: {key} must be nonempty text')
    reports = write_d_files(decisions, r_dir, d_dir, meta, args.dry_run, args.overwrite)
    for report in reports:
        print(f"{report['status']}: {report['output']}")
        if args.dry_run:
            for decision in report['decisions']:
                print(f"  profile {decision['target']['profile_index']}: PRES/TEMP/PSAL -> bad (4)")
                for suggestion in decision['suggestions']:
                    print(f"    {suggestion['checker']}: {suggestion['reason']}")
    return reports


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, KeyError) as exc:
        print(f'DMQC error: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
