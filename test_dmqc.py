"""Isolated tests for instruction merging and NetCDF output; no network needed."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset
import yaml

import dmqc_process as cli
from dmqc.download import download_r_files
from dmqc.combine import combine_instructions
from dmqc.instructions import (load_instructions,
                               parse_document, sha256)
from dmqc.writer import put_text, text, write_d_files


NAME = 'R6903708_001.nc'


def make_source(path, adjusted=True, fixed_history=False):
    with Dataset(path, 'w', format='NETCDF3_CLASSIC') as ds:
        for name, size in {'N_PROF': 2, 'N_LEVELS': 4, 'N_PARAM': 3, 'N_CALIB': 1,
                           'N_HISTORY': 1 if fixed_history else None,
                           'STRING4': 4, 'STRING8': 8, 'STRING16': 16,
                           'STRING64': 64, 'STRING256': 256, 'DATE_TIME': 14}.items():
            ds.createDimension(name, size)

        def char(name, dims, values=None):
            v = ds.createVariable(name, 'S1', dims, fill_value=b' ')
            if values is not None:
                v[:] = values
            return v

        platform = char('PLATFORM_NUMBER', ('N_PROF', 'STRING8'))
        for ip in range(2):
            put_text(platform, ip, '6903708')
        ds.createVariable('CYCLE_NUMBER', 'i4', ('N_PROF',))[:] = [1, 1]
        char('DATA_MODE', ('N_PROF',), np.array([b'A', b'A']))
        dsi = char('DATA_STATE_INDICATOR', ('N_PROF', 'STRING4'))
        for ip in range(2):
            put_text(dsi, ip, '2B')
        put_text(char('DATE_UPDATE', ('DATE_TIME',)), slice(None), '20200101000000')
        station = char('STATION_PARAMETERS', ('N_PROF', 'N_PARAM', 'STRING16'))
        char('PARAMETER_DATA_MODE', ('N_PROF', 'N_PARAM'), np.full((2, 3), b'A'))
        parameter = char('PARAMETER', ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING16'))
        for j, p in enumerate(('PRES', 'TEMP', 'PSAL')):
            for ip in range(2):
                put_text(station, (ip, j), p)
                put_text(parameter, (ip, 0, j), p)
            for suffix in ('', '_ADJUSTED', '_ADJUSTED_ERROR'):
                v = ds.createVariable(p + suffix, 'f4', ('N_PROF', 'N_LEVELS'), fill_value=99999.)
                values = np.ma.array([[1, 2, 3, 0], [5, 6, 7, 0]], mask=[[0, 0, 0, 1]] * 2)
                if suffix == '' or adjusted:
                    v[:] = values
            char(p + '_QC', ('N_PROF', 'N_LEVELS'), np.array([[b'1', b'3', b'4', b'9']] * 2))
            char(p + '_ADJUSTED_QC', ('N_PROF', 'N_LEVELS'),
                 np.array([[b'1', b'3', b'4', b'9']] * 2) if adjusted else None)
            char('PROFILE_' + p + '_QC', ('N_PROF',), np.array([b'A', b'A']))
        for field in ('EQUATION', 'COEFFICIENT', 'COMMENT', 'DATE'):
            v = char('SCIENTIFIC_CALIB_' + field,
                     ('N_PROF', 'N_CALIB', 'N_PARAM', 'DATE_TIME' if field == 'DATE' else 'STRING256'))
            for ip in range(2):
                for j in range(3):
                    put_text(v, (ip, 0, j), '20200101000000' if field == 'DATE' else 'original')
        for field, width in {'DATE': 'DATE_TIME', 'ACTION': 'STRING4', 'STEP': 'STRING4',
                             'PARAMETER': 'STRING16', 'INSTITUTION': 'STRING4',
                             'SOFTWARE': 'STRING4', 'SOFTWARE_RELEASE': 'STRING4',
                             'REFERENCE': 'STRING64'}.items():
            v = char('HISTORY_' + field, ('N_HISTORY', 'N_PROF', width))
            for ip in range(2):
                put_text(v, (0, ip), '20200101000000' if field == 'DATE' else 'old')


def document(action='flag', index=0, checker='test_checker'):
    entry = {'target': {'source': NAME, 'profile_index': index, 'selection': 'whole_profile'},
             'action': action, 'reason': 'Synthetic bad-profile check'}
    if action == 'flag':
        entry['flag'] = '4'
    return {'schema_version': 1, 'checker': checker, 'instructions': [entry]}


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.r_dir = self.root / '6903708' / 'R'
        self.r_dir.mkdir(parents=True)
        self.source = self.r_dir / NAME
        make_source(self.source)
        self.d_dir = self.r_dir.parent / 'D'

    def decisions(self, *docs):
        instructions = []
        for i, doc in enumerate(docs or (document(),)):
            instructions.extend(parse_document(doc, f'checker-{i}.yaml'))
        return combine_instructions(instructions, self.r_dir)

    def test_bad_wins_in_either_order_and_keeps_provenance(self):
        bad, good = document(), document('no_finding', checker='other')
        for docs in ((bad, good), (good, bad)):
            decisions = self.decisions(*docs)
            self.assertEqual(len(decisions), 1)
            self.assertEqual(len(decisions[0].suggestions), 2)
        self.assertEqual(self.decisions(good), [])

    def test_writer_masks_adjusted_values_preserves_raw_and_other_profile(self):
        original_hash = sha256(self.source)
        reports = write_d_files(self.decisions(), self.r_dir, self.d_dir)
        self.assertEqual(sha256(self.source), original_hash)
        output = self.d_dir / ('D' + NAME[1:])
        with Dataset(output) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'A'])
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED_QC'][0], [b'4', b'4', b'4', b'9'])
            self.assertEqual(ds['TEMP_ADJUSTED'][0].count(), 0)
            self.assertEqual(ds['TEMP_ADJUSTED_ERROR'][0].count(), 0)
            self.assertEqual(ds['TEMP_ADJUSTED'][1].count(), 3)
            self.assertEqual(ds['PROFILE_TEMP_QC'][0], b'F')
            self.assertEqual(len(ds.dimensions['N_CALIB']), 2)
            self.assertEqual(len(ds.dimensions['N_HISTORY']), 4)
            self.assertEqual(text(ds['HISTORY_STEP'][1, 0]), 'ARSQ')
            self.assertEqual(text(ds['SCIENTIFIC_CALIB_EQUATION'][0, 0, 1]), 'original')
            self.assertIn('test_checker', text(ds['SCIENTIFIC_CALIB_COMMENT'][0, 1, 1]))
            self.assertEqual(text(ds['DATE_UPDATE'][:])[:2], '20')
        audit = json.loads(output.with_suffix('.report.json').read_text())
        self.assertEqual(audit['source_sha256'], original_hash)
        self.assertEqual(audit['output_sha256'], sha256(output))
        self.assertEqual(reports[0]['changes'][0]['bad_samples'], 3)

    def test_raw_only_profile_and_nonzero_profile_target(self):
        make_source(self.source, adjusted=False)
        write_d_files(self.decisions(document(index=1)), self.r_dir, self.d_dir)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'A', b'D'])
            np.testing.assert_array_equal(ds['PSAL_ADJUSTED_QC'][1], [b'4', b'4', b'4', b'9'])

    def test_both_profiles_share_one_output_and_fixed_dimensions_grow(self):
        make_source(self.source, fixed_history=True)
        reports = write_d_files(self.decisions(document(), document(index=1)), self.r_dir, self.d_dir)
        self.assertEqual(len(reports), 1)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])
            self.assertEqual(len(ds.dimensions['N_HISTORY']), 7)

    def test_dry_run_and_no_finding_write_nothing(self):
        write_d_files(self.decisions(), self.r_dir, self.d_dir, dry_run=True)
        write_d_files(self.decisions(document('no_finding')), self.r_dir, self.d_dir)
        self.assertFalse(self.d_dir.exists())

    def test_stale_source_and_invalid_profile_are_rejected(self):
        decisions = self.decisions()
        with Dataset(self.source, 'r+') as ds:
            ds['TEMP'][0, 0] = 20
        with self.assertRaisesRegex(ValueError, 'source changed'):
            write_d_files(decisions, self.r_dir, self.d_dir)
        with self.assertRaisesRegex(ValueError, 'profile index'):
            write_d_files(self.decisions(document(index=7)), self.r_dir, self.d_dir)
        self.assertFalse(self.d_dir.exists())

    def test_explicit_checksum_mismatch_is_rejected(self):
        doc = document()
        doc['instructions'][0]['source_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'source changed'):
            self.decisions(doc)

    def test_unsupported_actions_flags_selections_and_pending_review(self):
        for key, value in [('action', 'adjust'), ('flag', '1'), ('status', 'pending')]:
            doc = document()
            doc['instructions'][0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.decisions(doc)
        for key, value in [('selection', 'pressure_range'), ('profile_index', -1),
                           ('source', '../' + NAME), ('parameters', ['TEMP'])]:
            doc = document()
            doc['instructions'][0]['target'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.decisions(doc)

    def test_existing_output_requires_overwrite(self):
        decisions = self.decisions()
        write_d_files(decisions, self.r_dir, self.d_dir)
        with self.assertRaises(FileExistsError):
            write_d_files(decisions, self.r_dir, self.d_dir)
        self.assertEqual(len(write_d_files(decisions, self.r_dir, self.d_dir, overwrite=True)), 1)

    def test_validation_failure_publishes_nothing(self):
        with patch('dmqc.writer.validate_output', side_effect=ValueError('failed check')):
            with self.assertRaisesRegex(ValueError, 'failed check'):
                write_d_files(self.decisions(), self.r_dir, self.d_dir)
        self.assertEqual(list(self.d_dir.iterdir()), [])

    def test_legacy_good_is_no_finding_bad_requires_reason(self):
        legacy = self.r_dir.parent / 'cycles'
        legacy.mkdir()
        path = legacy / '001.yaml'
        path.write_text("cycle: '001'\nqc_flag: '1'\nnote: ''\n")
        loaded = load_instructions(self.root / 'instructions', self.r_dir, '6903708', legacy_dir=legacy)
        self.assertEqual(combine_instructions(loaded, self.r_dir), [])
        path.write_text("cycle: '001'\nqc_flag: '4'\nnote: 'Bad cast'\n")
        loaded = load_instructions(self.root / 'instructions', self.r_dir, '6903708', legacy_dir=legacy)
        self.assertEqual(len(combine_instructions(loaded, self.r_dir)), 1)

    def test_cli_defaults_preserve_current_float_and_all_stage(self):
        args = cli.parser().parse_args([])
        self.assertEqual(args.float_id, '6903708')
        self.assertEqual(args.stage, 'all')
        self.assertEqual(args.work_dir, cli.DEFAULT_WORK_DIR)
        with patch('dmqc_process.download_r_files') as download, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(['--work-dir', str(self.root)]), 0)
            self.assertEqual(download.call_args.args[0], '6903708')

    def test_cli_write_from_partner_yaml_and_cycle_filter(self):
        folder = self.r_dir.parent / 'instructions' / 'partner'
        folder.mkdir(parents=True)
        (folder / 'check.yaml').write_text(yaml.safe_dump(document()))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(['write', '--work-dir', str(self.root), '--cycles', '2']), 0)
            self.assertFalse(self.d_dir.exists())
            self.assertEqual(cli.main(['write', '--work-dir', str(self.root)]), 0)
        self.assertTrue((self.d_dir / ('D' + NAME[1:])).is_file())

    def test_download_is_atomic_and_skips_existing_sources(self):
        payload = self.source.read_bytes()
        directory = self.root / 'download'
        def open_url(url, **kwargs):
            return io.BytesIO(f'<a href="{NAME}">file</a>'.encode() if url.endswith('/') else payload)
        with patch('dmqc.download.urlopen', side_effect=open_url):
            download_r_files('6903708', directory, dry_run=True)
            self.assertFalse(directory.exists())
            download_r_files('6903708', directory)
            checksum = sha256(directory / NAME)
            download_r_files('6903708', directory)
            self.assertEqual(sha256(directory / NAME), checksum)
        failed = self.root / 'failed-download'
        with patch('dmqc.download.urlopen', side_effect=[io.BytesIO(f'<a href="{NAME}">'.encode()), io.BytesIO(b'not netcdf')]):
            with self.assertRaises(OSError):
                download_r_files('6903708', failed)
        self.assertEqual(list(failed.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
