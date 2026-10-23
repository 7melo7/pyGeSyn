import copy
import unittest

from pyGeSyn.find import _find_candidates, _pick_best


def h(qs, qe, ss, se, strand='+', chrom='chr6'):
    return dict(qstart=qs,qend=qe,sstart=ss,send=se,strand=strand,chrom=chrom,
                length=qe-qs+1,identity=98,bitscore=1000)


class FlankChainTests(unittest.TestCase):
    def test_fragmented_short_terminal_support_and_outside_repeats(self):
        hits=[h(1,20138,10674654,10694784),
              h(126589,136103,10715425,10724948),
              h(136090,139062,10725035,10728007),
              h(138851,143179,10727887,10732228),
              h(143177,146528,10732270,10735646),
              h(163818,166082,10740358,10742625),
              h(166039,167306,10742623,10743880),
              h(175217,178704,10746213,10749691),
              h(72532,73819,10557287,10558577,'-'),
              h(108047,109438,10557287,10558680,'-')]
        before=copy.deepcopy(hits)
        best=_pick_best(_find_candidates(hits,178704,0))
        self.assertEqual((best['start'],best['end']),(10674654,10749691))
        self.assertGreater(best['right_chain_hits'],1)
        self.assertEqual(best['support_status'],'paired_flank_chains')
        self.assertEqual(hits,before)
        reverse=[dict(x,sstart=40000000-x['send'],send=40000000-x['sstart'],
                      strand='-' if x['strand']=='+' else '+') for x in hits]
        rb=_pick_best(_find_candidates(reverse,178704,0))
        self.assertEqual((rb['start'],rb['end']),(40000000-best['end'],40000000-best['start']))

    def test_ambiguous_copies_are_reported(self):
        hits=[h(1,15000,100001,115000),h(85001,100000,185001,200000)]
        hits += [dict(x,sstart=x['sstart']+300000,send=x['send']+300000) for x in hits]
        best=_pick_best(_find_candidates(hits,100000,0))
        self.assertEqual(best['support_status'],'ambiguous_pair')
        self.assertAlmostEqual(best['score_margin'],0)

    def test_fallback_and_weak_support_are_explicit(self):
        best=_pick_best(_find_candidates([h(30000,40000,1000,11000),h(50000,60000,15000,25000)],100000,0))
        self.assertEqual(best['support_status'],'legacy_unresolved')
        best=_pick_best(_find_candidates([h(1,100,1000,1099),h(90001,100000,50000,59999)],100000,0))
        self.assertEqual(best['support_status'],'weak_flank_support')

    def test_short_queries_and_alias(self):
        hits=[h(1,3,10,12)]
        for mode in ('chains','anchors'):
            best=_pick_best(_find_candidates(hits,3,0,min_hits=1,boundary_mode=mode))
            self.assertEqual((best['start'],best['end']),(10,12))


if __name__=='__main__':
    unittest.main()
