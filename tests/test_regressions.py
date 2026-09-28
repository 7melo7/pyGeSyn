import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from pyGeSyn.cli import _init_workdir, main, parse_regions_file
from pyGeSyn.dotplot import _blast_pair
from pyGeSyn.find import _find_candidates, discover_regions, write_results


def hit(qs, qe, ss, se, strand='+', identity=99):
    return dict(chrom='chr1', qstart=qs, qend=qe, sstart=ss, send=se,
                strand=strand, identity=identity, length=qe-qs+1, bitscore=800)


class RegressionTests(unittest.TestCase):
    def test_reverse_candidates_roundtrip(self):
        hsps = [hit(1, 400, 1601, 2000, '-'), hit(601, 1000, 1001, 1400, '-')]
        c = _find_candidates(hsps, 1000, 0)[0]
        self.assertEqual((c['start'], c['end'], c['subject_span']), (1001, 2000, 1000))
        self.assertEqual(c['orientation'], '-')
        self.assertEqual(c['identity'], 99)
        c['genome'] = 'B'
        with tempfile.TemporaryDirectory() as tmp:
            best = Path(tmp) / 'best.csv'
            write_results([c], [c], Path(tmp) / 'all.tsv', best)
            self.assertEqual(parse_regions_file(best), [('B', 'chr1', 1001, 2000)])

    def test_mixed_orientation_and_weighted_identity(self):
        c = _find_candidates([hit(1, 100, 1, 100, identity=80),
                              hit(101, 400, 101, 400, '-', 100)], 400, 0)[0]
        self.assertEqual(c['orientation'], 'mixed')
        self.assertEqual(c['identity'], 95)

    def test_min_hits(self):
        hsps = [hit(1, 400, 1, 400), hit(601, 1000, 601, 1000)]
        self.assertFalse(_find_candidates(hsps, 1000, 0, min_hits=3))
        self.assertTrue(_find_candidates(hsps[:1], 1000, 0, min_hits=1))
        self.assertFalse(_find_candidates(hsps[:1], 1000, 0))

    def test_invalid_parameters(self):
        for options in ({'min_hits': 0}, {'window_mult': 0},
                        {'window_mult': float('nan')}, {'min_coverage': 2}):
            kwargs = dict(query_len=1000, min_coverage=0)
            kwargs.update(options)
            with self.assertRaises(ValueError):
                _find_candidates([], **kwargs)
        self.assertTrue(_find_candidates([hit(1, 1, 1, 1)], 1, 0,
                                         window_mult=0.1, min_hits=1))

    def test_removed_gap_options_fail_explicitly(self):
        for key in ('cluster_gap', 'merge_gap'):
            with self.assertRaisesRegex(ValueError, 'unsupported'):
                discover_regions('unused', 'unused', **{key: 100})

    def test_serial_parallel_min_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / 'a.fa'
            fasta.write_text('>chr1\n' + 'A' * 1000)
            cfg = Path(tmp) / 'cfg.json'
            cfg.write_text(json.dumps({n: {'fasta': str(fasta)} for n in ('A', 'B', 'C')}))
            hsps = [hit(1, 400, 1, 400), hit(601, 1000, 601, 1000)]
            with patch('pyGeSyn.find._blast_genome', return_value=hsps):
                for threads in (1, 2):
                    for threshold, expected in ((2, 2), (3, 0)):
                        _, best = discover_regions('A$chr1$1$1000', str(cfg),
                                                   min_hits=threshold, threads=threads)
                        self.assertEqual(len(best), expected)

    def test_workdir_preserves_existing_files(self):
        previous = Path.cwd()
        with tempfile.TemporaryDirectory() as tmp:
            try:
                os.chdir(tmp)
                Path('temp/nested').mkdir(parents=True)
                sentinel = Path('temp/keep.txt')
                sentinel.write_text('keep')
                first, second = _init_workdir(), _init_workdir()
                self.assertNotEqual(first, second)
                self.assertTrue(first.is_dir() and second.is_dir())
                self.assertEqual(sentinel.read_text(), 'keep')
                self.assertTrue(Path('temp/nested').is_dir())
            finally:
                os.chdir(previous)

    def test_cli_defaults(self):
        with patch('sys.argv', ['pyGeSyn', 'find', 'query.fa', 'config.json']), \
             patch('pyGeSyn.cli.cmd_find') as run:
            main()
            args = run.call_args.args[0]
            self.assertEqual((args.min_length, args.min_identity), (200, 80))

    def test_dotplot_retains_reverse_coordinates(self):
        result = subprocess.CompletedProcess([], 0, '1\t2000\t2000\t1\t2000\t100\n', '')
        with patch('pyGeSyn.dotplot.subprocess.run', return_value=result):
            self.assertEqual(_blast_pair('ACGT', 'ACGT', 100, 80)[0]['send'], 1)

    def test_blast_errors_are_not_empty_results(self):
        result = subprocess.CompletedProcess([], 1, '', 'invalid query')
        with patch('pyGeSyn.dotplot.subprocess.run', return_value=result):
            with self.assertRaisesRegex(RuntimeError, 'invalid query'):
                _blast_pair('ACGT', 'ACGT', 100, 80)

    @unittest.skipUnless(shutil.which('blastn'), 'BLAST+ is not installed')
    def test_real_reverse_complement(self):
        rng = random.Random(42)
        seq = ''.join(rng.choices('ACGT', k=2000))
        reverse = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        hsps = _blast_pair(seq, reverse, 100, 80)
        self.assertTrue(hsps)
        self.assertTrue(all((h['qend']-h['qstart'])*(h['send']-h['sstart']) < 0 for h in hsps))

    @unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),
                         'BLAST+ is not installed')
    def test_cli_reverse_pipeline(self):
        rng = random.Random(7)
        seq = ''.join(rng.choices('ACGT', k=2000))
        reverse = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = {}
            for name, sequence in (('A', seq), ('B', reverse)):
                fa = root / f'{name}.fa'
                fa.write_text(f'>chr1\n{sequence}\n')
                config[name] = {'fasta': str(fa)}
            (root / 'genomes.json').write_text(json.dumps(config))
            (root / 'temp/nested').mkdir(parents=True)
            sentinel = root / 'temp/keep.txt'
            sentinel.write_text('keep')
            env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
            commands = [
                ['find', 'A$chr1$1$2000', 'genomes.json', '--min-hits', '1'],
                ['plot', 'genomes.json', 'best_regions.csv', '-o', 'plot.png'],
                ['dotplot', 'best_regions.csv', 'genomes.json', 'A,B', '-o', 'dotplot.png'],
                ['dedup', 'best_regions.csv', 'genomes.json'],
            ]
            for command in commands:
                result = subprocess.run([sys.executable, '-B', '-m', 'pyGeSyn', *command],
                                        cwd=tmp, env=env, capture_output=True, text=True,
                                        timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(parse_regions_file(root / 'best_regions.csv'),
                             [('A', 'chr1', 1, 2000), ('B', 'chr1', 1, 2000)])
            self.assertEqual(len(parse_regions_file(root / 'nonredundant.csv')), 1)
            for output in ('plot.png', 'dotplot.png'):
                self.assertGreater((root / output).stat().st_size, 1000)
            self.assertEqual(sentinel.read_text(), 'keep')
            self.assertEqual(len(list((root / 'temp').glob('pyGeSyn-*'))), 3)


if __name__ == '__main__':
    unittest.main()
