"""Inspect core Argo D-files and save a YAML inventory and validation report.

Accepts a float directory or its D/ directory. With no arguments, checks 6903708.
Files are opened read-only. This implements a documented subset of Argo rules;
it does not replace the official GDAC format checker or scientific review.
"""
import argparse
import os
from pathlib import Path
import sys
import tempfile

import yaml

from dmqc.verification import verify_directory


DEFAULT_DIRECTORY = Path(r'C:\Data\ARGO_Dataa\DMQCprocessing\6903708' if os.name == 'nt'
                         else '/mnt/c/Data/ARGO_Dataa/DMQCprocessing/6903708')


def save_report(report, output):
    """Replace the YAML report atomically, leaving an earlier report on failure."""
    output = Path(output)
    if output.suffix.lower() not in ('.yaml', '.yml'):
        raise ValueError('Report output must have a .yaml or .yml extension')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=output.parent,
                                         suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            yaml.safe_dump(report, stream, sort_keys=False, allow_unicode=True)
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=DEFAULT_DIRECTORY,
                        help='Float folder or D folder (default: %(default)s)')
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--output', type=Path, help='Report path; default: reports/verification.yaml')
    output.add_argument('--no-save', action='store_true', help='Print results without writing a report')
    parser.add_argument('--verbose', action='store_true', help='Print each finding as well as the file summary')
    parser.add_argument('--strict', action='store_true', help='Return failure for warnings too')
    args = parser.parse_args(argv)
    try:
        report = verify_directory(args.directory)
        for result in report['files']:
            parameters = sorted({p for profile in result['profiles'] for p in profile.get('station_parameters', [])})
            print(f"{result['status']:4} {result['file']}: {len(result['profiles'])} profiles, "
                  f"{', '.join(parameters) or 'unknown parameters'}; "
                  f"{result['error_count']} errors, {result['warning_count']} warnings")
            if args.verbose:
                for issue in result['issues']:
                    location = issue.get('variable', '')
                    if 'profile_index' in issue:
                        location += f" profile {issue['profile_index']}"
                    count = f" ({issue['count']} samples)" if 'count' in issue else ''
                    print(f"  {issue['severity']} {issue['code']} {location}: {issue['message']}{count}")
        summary = report['summary']
        print(f"\n{summary['files']} files / {summary['profiles']} profiles: "
              f"{summary['passed']} passed, {summary['failed']} failed, "
              f"{summary['warning_only']} with warnings only.")
        for code, count in summary['issue_counts'].items():
            print(f'  {code}: {count} findings')
        for issue in report['directory_issues']:
            print(f"  {issue['message']}: {', '.join(issue['files'])}")
        print(report['scope'])
        if not args.no_save:
            directory = Path(report['directory'])
            root = directory.parent if directory.name == 'D' else directory
            destination = args.output or root / 'reports' / 'verification.yaml'
            save_report(report, destination)
            print(f'Report: {destination}')
        return int(bool(summary['errors'] or (args.strict and summary['warnings'])))
    except (ValueError, OSError, yaml.YAMLError) as exc:
        print(f'Verification error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
