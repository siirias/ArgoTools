"""Check regional surface salinity and write whole-profile rejection instructions."""
import argparse
from pathlib import Path
import sys

from dmqc.settings import add_settings_argument, float_directory
from dmqc.surface_salinity import load_config, check_directory, save_instructions
from verify_dfiles import save_report

DEFAULT_CONFIG = Path(__file__).resolve().parent / 'config' / 'surface_salinity.yaml'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_settings_argument(parser)
    parser.add_argument('directory', nargs='?', type=Path,
                        help='Float folder or R folder (default: local settings)')
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--dry-run', action='store_true', help='Report counts without saving files')
    args = parser.parse_args(argv)
    try:
        args.directory = float_directory(args.directory, args.settings)
        config = load_config(args.config)
        r_dir, entries, report = check_directory(args.directory, config)
        report['configuration_path'] = str(args.config.resolve())
        root = r_dir.parent if r_dir.name == 'R' else r_dir
        summary = report['summary']
        print(f"{len(report['sources'])} files, {len(report['profiles'])} profiles: "
              f"{summary.get('pass', 0)} passed, {summary.get('reject', 0)} rejected, "
              f"{summary.get('not_evaluated', 0)} not evaluated")
        for reason, count in report['skip_reasons'].items():
            print(f'  {reason}: {count}')
        if not args.dry_run:
            instructions = root / 'instructions' / 'surface_salinity.yaml'
            output = root / 'reports' / 'surface_salinity.yaml'
            save_report(report, output)
            save_instructions(instructions, entries)
            print(f'Instructions: {instructions}\nReport: {output}')
        return 0
    except (ValueError, OSError, KeyError, RuntimeError) as exc:
        print(f'Surface salinity check failed: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
