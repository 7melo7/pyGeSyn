import contextlib
import io
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

from pyGeSyn.cli import main, parse_regions_file
from pyGeSyn.find import (_apply_repeatmask, _load_te_index, _query_te_context,
                          _strictly_in_te, discover_regions)


def hsp(qs, qe, ss, se, chrom='chr1'):
    return dict(qstart=qs, qend=qe, sstart=ss, send=se, chrom=chrom,
                strand='+', identity=99, length=qe-qs+1, bitscore=100)


class RepeatMaskTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / 'genome.fa'
        self.fasta.write_text('>chr1\n' + 'A' * 10000)
        self.query_te = self.te('query.gff', [(1100, 1201)])
        self.subject_te = self.te('subject.gff', [(1000, 2000)])

    def te(self, name, regions):
        path = self.root / name
        path.write_text(''.join(f'chr1\t.\tTE\t{s}\t{e}\t.\t+\t.\tID=t{i}\n'
                                for i, (s, e) in enumerate(regions)))
        return str(path)

    def test_strict_boundaries_reverse_and_chromosome(self):
        index = _load_te_index(self.subject_te)
        for start, end, expected in ((1001, 1999, True), (1999, 1001, True),
                                     (1000, 1999, False), (1001, 2000, False),
                                     (999, 2001, False), (900, 1100, False)):
            self.assertEqual(_strictly_in_te(index['chr1'], start, end), expected)
        self.assertFalse(_strictly_in_te(index.get('chr2'), 1100, 1200))

    def test_single_te_not_union_and_nested_intervals(self):
        path = self.te('multiple.gff', [(1, 100), (90, 200), (30, 40)])
        intervals = _load_te_index(path)['chr1']
        self.assertFalse(_strictly_in_te(intervals, 80, 120))
        self.assertTrue(_strictly_in_te(intervals, 50, 80))

    def test_query_offset_or_subject_masking_and_mode_zero(self):
        cfg = {'A': {'te': self.query_te}}
        query_te, offset = _query_te_context('A$chr1$1001$3000', cfg)
        self.assertEqual(offset, 1000)
        hits = [hsp(101, 200, 5000, 5100), hsp(301, 400, 1101, 1200),
                hsp(501, 600, 1000, 1500), hsp(601, 700, 1500, 2000)]
        paths = {'fasta': str(self.fasta), 'te': self.subject_te}
        self.assertEqual(_apply_repeatmask(hits, paths, 1, query_te, offset), hits[2:])
        with patch('pyGeSyn.find._load_te_index', side_effect=AssertionError('must not load')):
            self.assertEqual(_apply_repeatmask(hits, paths, 0, query_te, offset), hits)

    def test_fasta_and_missing_annotations(self):
        with contextlib.redirect_stderr(io.StringIO()):
            query_te, offset = _query_te_context(str(self.fasta), {})
            self.assertIsNone(query_te)
            self.assertEqual(_query_te_context('A$chr1$1$2000', {'A': {}}), (None, 0))
            hit = hsp(1, 100, 1100, 1200)
            self.assertFalse(_apply_repeatmask([hit], {'te': self.subject_te}, 1, query_te, offset))
            self.assertEqual(_apply_repeatmask([hit], {'fasta': str(self.fasta)}, 1, None, 0), [hit])

    def test_serial_parallel_modes(self):
        cfg = self.root / 'genomes.json'
        cfg.write_text(json.dumps({
            'A': {'fasta': str(self.fasta), 'te': self.query_te},
            'B': {'fasta': str(self.fasta), 'te': self.subject_te},
            'C': {'fasta': str(self.fasta), 'te': self.subject_te}}))
        hits = [hsp(101, 200, 5000, 5100), hsp(301, 400, 1101, 1200),
                hsp(501, 600, 1000, 1500)]
        with patch('pyGeSyn.find._blast_genome', return_value=hits):
            for threads in (1, 2):
                for mode in (0, 1):
                    with patch('pyGeSyn.find._find_candidates', return_value=[]) as find:
                        discover_regions('A$chr1$1001$3000', str(cfg),
                                         threads=threads, repeatmaskmode=mode)
                        self.assertEqual(find.call_count, 2)
                        for call in find.call_args_list:
                            self.assertEqual(call.args[0], hits if mode == 0 else hits[2:])
        with self.assertRaisesRegex(ValueError, 'repeatmaskmode'):
            discover_regions('unused', 'unused', repeatmaskmode=2)

    def test_cli_choices_and_default(self):
        for option, expected in (([], 0), (['--repeatmaskmode', '1'], 1)):
            with patch('sys.argv', ['pyGeSyn', 'find', 'q.fa', 'g.json', *option]), \
                 patch('pyGeSyn.cli.cmd_find') as run:
                main()
                self.assertEqual(run.call_args.args[0].repeatmaskmode, expected)
        with patch('sys.argv', ['pyGeSyn', 'find', 'q.fa', 'g.json', '--repeatmaskmode', '2']), \
             contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                main()

    @unittest.skipUnless(shutil.which('blastn'), 'BLAST+ is not installed')
    def test_real_blast_cli(self):
        rng = random.Random(14)
        sequence = ''.join(rng.choices('ACGT', k=2000))
        query = self.root / 'query.fa'
        query.write_text('>query\n' + sequence)
        self.fasta.write_text('>chr1\n' + 'N' * 100 + sequence + 'N' * 100)
        te = self.te('real.gff', [(100, 2200)])
        cfg = self.root / 'real.json'
        cfg.write_text(json.dumps({'B': {'fasta': str(self.fasta), 'te': te}}))
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
        for mode in (0, 1):
            best = self.root / f'best{mode}.csv'
            result = subprocess.run([sys.executable, '-B', '-m', 'pyGeSyn', 'find',
                                     str(query), str(cfg), '--repeatmaskmode', str(mode),
                                     '--min-hits', '1', '--best-output', str(best)],
                                    cwd=self.root, env=env, capture_output=True, text=True,
                                    timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(parse_regions_file(best)), 1 if mode == 0 else 0)


if __name__ == '__main__':
    unittest.main()
