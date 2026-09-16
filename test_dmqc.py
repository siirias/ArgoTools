"""Isolated tests for instruction merging and NetCDF output; no network needed."""
import contextlib
import io
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
                               parse_document, sha256, write_instructions, write_profile_flags, read_yaml)
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
            self.assertEqual(len(decisions), 2)
            self.assertEqual(sum(d.rejected for d in decisions), 1)
            self.assertEqual(len(decisions[0].suggestions), 2)
        self.assertFalse(any(d.rejected for d in self.decisions(good)))

    def test_writer_masks_rejections_and_preserves_raw_and_other_profile_qc(self):
        original_hash = sha256(self.source)
        reports = write_d_files(self.decisions(), self.r_dir, self.d_dir)
        self.assertEqual(sha256(self.source), original_hash)
        output = self.d_dir / ('D' + NAME[1:])
        with Dataset(output) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED_QC'][0], [b'4', b'4', b'4', b'9'])
            self.assertEqual(ds['TEMP_ADJUSTED'][0].count(), 0)
            self.assertEqual(ds['TEMP_ADJUSTED_ERROR'][0].count(), 0)
            self.assertEqual(ds['TEMP_ADJUSTED'][1].count(), 2)
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED_QC'][1], [b'1', b'3', b'4', b'9'])
            self.assertEqual(ds['PROFILE_TEMP_QC'][0], b'F')
            self.assertEqual(len(ds.dimensions['N_CALIB']), 2)
            self.assertEqual(len(ds.dimensions['N_HISTORY']), 7)
            self.assertEqual(text(ds['HISTORY_STEP'][1, 0]), 'ARSQ')
            self.assertEqual(text(ds['SCIENTIFIC_CALIB_EQUATION'][0, 0, 1]), 'original')
            self.assertIn('test_checker', text(ds['SCIENTIFIC_CALIB_COMMENT'][0, 1, 1]))
            self.assertEqual(text(ds['DATE_UPDATE'][:])[:2], '20')
        report_path = self.d_dir.parent / 'reports' / output.with_suffix('.report.yaml').name
        audit = yaml.safe_load(report_path.read_text())
        self.assertEqual(reports[0]['report'], str(report_path))
        self.assertEqual(list(self.d_dir.glob('*.report.*')), [])
        with Dataset(output) as ds:
            self.assertEqual(text(ds['HISTORY_REFERENCE'][1, 0]), '../reports/' + report_path.name)
        self.assertEqual(audit['source_sha256'], original_hash)
        self.assertEqual(audit['output_sha256'], sha256(output))
        self.assertEqual(reports[0]['changes'][0]['bad_samples'], 3)

    def test_raw_only_profile_and_nonzero_profile_target(self):
        make_source(self.source, adjusted=False)
        write_d_files(self.decisions(document(index=1)), self.r_dir, self.d_dir)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])
            np.testing.assert_array_equal(ds['PSAL_ADJUSTED_QC'][1], [b'4', b'4', b'4', b'9'])

    def test_both_profiles_share_one_output_and_fixed_dimensions_grow(self):
        make_source(self.source, fixed_history=True)
        reports = write_d_files(self.decisions(document(), document(index=1)), self.r_dir, self.d_dir)
        self.assertEqual(len(reports), 1)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])
            self.assertEqual(len(ds.dimensions['N_HISTORY']), 7)

    def test_dry_run_writes_nothing_and_no_finding_still_produces_output(self):
        write_d_files(self.decisions(), self.r_dir, self.d_dir, dry_run=True)
        self.assertFalse(self.d_dir.exists())
        self.assertFalse((self.d_dir.parent / 'reports').exists())
        reports = write_d_files(self.decisions(document('no_finding')), self.r_dir, self.d_dir)
        self.assertEqual(len(reports), 1)
        self.assertTrue(all(c['outcome'] == 'retain' for c in reports[0]['changes']))

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

    def test_cycle_directory_does_not_contribute_instructions(self):
        cycles = self.r_dir.parent / 'cycles'
        cycles.mkdir()
        (cycles / '001.yaml').write_text("cycle: '001'\nqc_flag: '4'\nnote: 'Old decision'\n")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(['write', '--work-dir', str(self.root)]), 0)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED_QC'][0], [b'1', b'3', b'4', b'9'])

    def test_combine_stdout_is_readable_yaml(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(cli.main(['combine', '--work-dir', str(self.root)]), 0)
        plan = yaml.safe_load(stream.getvalue())
        self.assertEqual(plan['float_id'], '6903708')
        self.assertEqual(plan['file_count'], 1)
        self.assertEqual(plan['profile_count'], 2)
        self.assertEqual(plan['decisions'][0]['target']['source'], NAME)

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

    def test_all_files_are_written_and_good_values_and_errors_are_retained(self):
        other = self.r_dir / 'R6903708_002.nc'
        make_source(other)
        with Dataset(other, 'r+') as ds:
            ds['CYCLE_NUMBER'][:] = [2, 2]
            for p in ('PRES', 'TEMP', 'PSAL'):
                ds[p + '_ADJUSTED_QC'][:] = np.array([[b'1', b'1', b'1', b'9']] * 2)
                ds[p + '_ADJUSTED_ERROR'][:, :3] = 0.02
        decisions = self.decisions()
        self.assertEqual(len(decisions), 4)
        reports = write_d_files(decisions, self.r_dir, self.d_dir)
        self.assertEqual(len(reports), 2)
        with Dataset(other) as src, Dataset(self.d_dir / 'D6903708_002.nc') as ds:
            for p in ('PRES', 'TEMP', 'PSAL'):
                for suffix in ('_ADJUSTED', '_ADJUSTED_QC', '_ADJUSTED_ERROR'):
                    np.testing.assert_array_equal(ds[p + suffix][:], src[p + suffix][:])
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])
            self.assertEqual(text(ds['SCIENTIFIC_CALIB_EQUATION'][0, 1, 0]), 'original')

    def test_empty_rejection_list_exports_raw_only_good_data_without_inventing_errors(self):
        make_source(self.source, adjusted=False)
        decisions = combine_instructions([], self.r_dir)
        self.assertEqual(len(decisions), 2)
        reports = write_d_files(decisions, self.r_dir, self.d_dir)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED_QC'][0], [b'1', b'3', b'4', b'9'])
            np.testing.assert_array_equal(ds['TEMP_ADJUSTED'][0, :2], [1, 2])
            self.assertTrue(ds['TEMP_ADJUSTED'][0, 2:].mask.all())
            self.assertEqual(ds['TEMP_ADJUSTED_ERROR'][:].count(), 0)
        self.assertGreater(sum(c['samples_without_uncertainty'] for c in reports[0]['changes']), 0)

    def test_removing_rejection_restores_source_data_on_overwrite(self):
        write_d_files(self.decisions(), self.r_dir, self.d_dir)
        write_d_files(combine_instructions([], self.r_dir), self.r_dir, self.d_dir, overwrite=True)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(float(ds['TEMP_ADJUSTED'][0, 0]), 1.)
            self.assertEqual(ds['TEMP_ADJUSTED_QC'][0, 0], b'1')

    def test_minimal_yaml_without_optional_fields(self):
        doc = document()
        del doc['instructions'][0]['reason']
        output = self.root / 'minimal.yaml'
        write_instructions(output, doc['checker'], doc['instructions'])
        loaded = parse_document(read_yaml(output), output)
        self.assertEqual(loaded[0].reason, '')
        self.assertIsNone(loaded[0].source_sha256)
        self.assertEqual(len(combine_instructions(loaded, self.r_dir)), 2)

    def test_saved_profile_flags_are_consumable_by_dfile_writer(self):
        folder = self.root / 'instructions'
        output = folder / 'visual_inspector.yaml'
        hashes = {NAME: sha256(self.source)}
        write_profile_flags(output, {(NAME, 0), (NAME, 1)}, source_hashes=hashes,
                            metadata={'operator': 'test'})
        loaded = load_instructions(folder, self.r_dir, '6903708')
        self.assertEqual(loaded[0].metadata, {'operator': 'test'})
        reports = write_d_files(combine_instructions(loaded, self.r_dir), self.r_dir, self.d_dir)
        self.assertEqual(len(reports), 1)
        with Dataset(self.d_dir / ('D' + NAME[1:])) as ds:
            self.assertEqual(ds['DATA_MODE'][:].tolist(), [b'D', b'D'])

    def test_yaml_replacement_is_atomic_and_empty_clears_own_report(self):
        folder = self.root / 'instructions'
        output = folder / 'visual_inspector.yaml'
        write_profile_flags(output, {(NAME, 0)})
        previous = output.read_bytes()
        partner = folder / 'partner.yaml'
        partner.write_bytes(previous)
        with patch('dmqc.instructions.os.replace', side_effect=OSError('disk error')):
            with self.assertRaises(OSError):
                write_profile_flags(output, set())
        self.assertEqual(output.read_bytes(), previous)
        self.assertEqual(sorted(p.name for p in folder.iterdir()), ['partner.yaml', 'visual_inspector.yaml'])
        write_profile_flags(output, set())
        self.assertEqual(read_yaml(output)['instructions'], [])
        self.assertEqual(partner.read_bytes(), previous)

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
