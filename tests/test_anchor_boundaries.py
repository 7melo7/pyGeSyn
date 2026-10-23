import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pyGeSyn.find import _find_candidates, _pick_best, discover_regions, write_results


def hit(qs, qe, ss, se, strand='+', chrom='chr1'):
    return dict(qstart=qs, qend=qe, sstart=ss, send=se, strand=strand,
                chrom=chrom, identity=99, length=qe-qs+1, bitscore=500)


def locus(reverse=False):
    if reverse:
        return [hit(1,1500,18501,20000,'-'), hit(1501,8500,11501,18500,'-'),
                hit(8501,10000,10001,11500,'-')]
    return [hit(1,1500,10001,11500), hit(1501,8500,11501,18500),
            hit(8501,10000,18501,20000)]


class AnchorTests(unittest.TestCase):
    def test_external_internal_repeats_do_not_change_bounds_or_score(self):
        for reverse in (False, True):
            baseline = _pick_best(_find_candidates(locus(reverse), 10000, 0))
            extra = [hit(3000,5000,21000,23000), hit(5000,7000,6000,8000)]
            best = _pick_best(_find_candidates(locus(reverse)+extra, 10000, 0))
            self.assertEqual((best['start'], best['end']), (10001,20000))
            self.assertEqual(best['boundary_method'], 'paired_flank_chains')
            self.assertEqual(best['anchor_strand'], '-' if reverse else '+')
            self.assertAlmostEqual(best['score'], baseline['score'])
            self.assertEqual(best['excluded_hsps'], 2)
            legacy = _pick_best(_find_candidates(locus(reverse)+extra, 10000, 0,
                                                boundary_mode='window'))
            self.assertEqual((legacy['start'], legacy['end']), (6000,23000))

    def test_local_inversion_is_retained(self):
        hits = locus()
        hits[1]['strand'] = '-'
        best = _pick_best(_find_candidates(hits, 10000, 0))
        self.assertEqual(best['orientation'], 'mixed')
        self.assertEqual(best['anchor_strand'], '+')
        self.assertEqual(best['coverage'], 1)

    def test_missing_or_weak_end_falls_back(self):
        for hits in (locus()[:2],):
            best = _pick_best(_find_candidates(hits,10000,0))
            self.assertEqual(best['boundary_method'],'window_fallback')

    def test_incompatible_end_strands_do_not_form_pair(self):
        hits = locus()
        hits[-1]['strand'] = '-'
        best = _pick_best(_find_candidates(hits,10000,0))
        self.assertEqual(best['boundary_method'],'window_fallback')

    def test_two_copies_not_selected_as_one_spanning_locus(self):
        hits = locus()
        second = [dict(h,sstart=h['sstart']+30000,send=h['send']+30000) for h in hits]
        best = _pick_best(_find_candidates(hits+second,10000,0))
        self.assertIn((best['start'],best['end']), [(10001,20000),(40001,50000)])

    def test_real_expansion_allowed(self):
        hits = [hit(1,1500,10001,11500),hit(8501,10000,38501,40000)]
        best = _pick_best(_find_candidates(hits,10000,0))
        self.assertEqual(best['boundary_method'],'paired_flank_chains')
        self.assertEqual(best['subject_span'],30000)
        self.assertGreater(best['span_penalty'],0)

    def test_single_full_length_hsp_and_min_hits(self):
        hits = [hit(1,10000,10001,20000)]
        best = _pick_best(_find_candidates(hits,10000,0,min_hits=1))
        self.assertEqual(best['boundary_method'],'paired_flank_chains')
        self.assertFalse(_find_candidates(hits,10000,0,min_hits=2))

    def test_serial_parallel_and_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fa=root/'a.fa'
            fa.write_text('>chr1\n'+'A'*60000)
            cfg=root/'genomes.json'
            cfg.write_text(json.dumps({n:{'fasta':str(fa)} for n in ('A','B','C')}))
            hits=locus()+[hit(3000,5000,21000,23000)]
            with patch('pyGeSyn.find._blast_genome',return_value=hits):
                for threads in (1,2):
                    for mode in ('chains','anchors','window'):
                        all_hits,best=discover_regions('A$chr1$1$10000',str(cfg),
                                                       threads=threads,boundary_mode=mode)
                        self.assertEqual(len(best),2)
                        self.assertTrue(all(c['end']==(20000 if mode!='window' else 23000) for c in best))
            write_results(all_hits,best,root/'all.tsv',root/'best.csv')
            with open(root/'all.tsv') as stream:
                rows=list(csv.DictReader(stream,delimiter='\t'))
            self.assertTrue(rows)
            self.assertEqual(rows[0]['boundary_method'],'window')
            self.assertIn('left_anchor_sstart',rows[0])


if __name__=='__main__':
    unittest.main()
