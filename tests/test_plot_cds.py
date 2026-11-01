import csv
import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pyGeSyn.plot_cds import representatives, compare, best_chain, reverse_complement


class PlotCdsTests(unittest.TestCase):
    def test_representative_full_sequence_negative_and_missing_exons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fa=root/'a.fa';fa.write_text('>c\n'+'AACCGGTT'*100)
            gff=root/'a.gff'
            gff.write_text('c\t.\tgene\t1\t100\t.\t-\t.\tID=g\n'
                          'c\t.\tmRNA\t1\t100\t.\t-\t.\tID=short;Parent=g\n'
                          'c\t.\tCDS\t1\t5\t.\t-\t0\tParent=short\n'
                          'c\t.\tmRNA\t1\t100\t.\t-\t.\tID=long;Parent=g\n'
                          'c\t.\tCDS\t1\t10\t.\t-\t2\tParent=long\n'
                          'c\t.\tCDS\t91\t100\t.\t-\t0\tParent=long\n'
                          'c\t.\texon\t1\t10\t.\t-\t.\tParent=long\n')
            colors=root/'colors';colors.write_text('short #FF0000')
            selected,features,manifest,_=representatives(dict(fasta=str(fa),gff3=str(gff),gene_colors=str(colors)),
                                                         ('A','c',5,95),root,0)
            self.assertEqual(selected[0]['transcript_id'],'long')
            seq='AACCGGTT'*100
            self.assertEqual(selected[0]['seq'],reverse_complement(seq[:10]+seq[90:100]))
            self.assertEqual(features[0]['color'],'#FF0000')
            self.assertEqual(len(features[0]['exons']),2)
            self.assertTrue(manifest[0]['boundary_clipped'])

    def test_overlap_and_reverse_chains(self):
        def h(qs,qe,ss,se,bits,strand='+'):
            return dict(qs=qs,qe=qe,ss=ss,se=se,bits=bits,strand=strand)
        a=h(1,50,1,50,100);b=h(40,90,40,90,80);c=h(51,100,51,100,90)
        self.assertEqual(best_chain([a,b,c]),[a,c])
        a=h(1,50,51,100,100,'-');c=h(51,100,1,50,90,'-')
        self.assertEqual(best_chain([c,a]),[a,c])

    def test_either_side_coverage_and_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            a=[dict(id='a',seq='A'*100)];b=[dict(id='b',seq='A'*1000)]
            result=subprocess.CompletedProcess([],0,'a\tb\t1\t90\t101\t190\t90\t81\t100\n','')
            with patch('pyGeSyn.plot_cds.subprocess.run',return_value=result):
                match=compare(a,b,root/'a',root/'b',root,0,90,.9)[0]
                self.assertTrue(match['passed'])
                self.assertEqual(match['subject_coverage'],.09)
                self.assertFalse(compare(a,b,root/'a',root/'b',root,0,91,.9)[0]['passed'])
                self.assertFalse(compare(a,b,root/'a',root/'b',root,0,90,.91)[0]['passed'])
            result.stdout='a\tb\t101\t190\t1\t90\t90\t81\t100\n'
            a[0]['seq']='A'*1000; b[0]['seq']='A'*100
            with patch('pyGeSyn.plot_cds.subprocess.run',return_value=result):
                match=compare(a,b,root/'a',root/'b',root,0,90,.9)[0]
                self.assertTrue(match['passed'])
                self.assertEqual(match['query_coverage'],.09)

    @unittest.skipUnless(shutil.which('blastn'),'Requires BLAST')
    def test_real_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq=''.join(random.Random(5).choices('ACGT',k=600))
            config={}
            for name in ('A','B'):
                fa=root/f'{name}.fa';fa.write_text('>c\n'+seq)
                gff=root/f'{name}.gff';gff.write_text('c\t.\tgene\t1\t600\t.\t+\t.\tID=g\n'
                    'c\t.\tmRNA\t1\t600\t.\t+\t.\tID=t;Parent=g\n'
                    'c\t.\tCDS\t1\t600\t.\t+\t0\tParent=t\n')
                config[name]=dict(fasta=str(fa),gff3=str(gff))
            cfg=root/'cfg.json';cfg.write_text(json.dumps(config));regions=root/'regions.csv'
            regions.write_text('A,c,1,600\nB,c,1,600\n')
            import os
            env=dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]))
            result=subprocess.run([sys.executable,'-B','-m','pyGeSyn','plot_cds',str(cfg),str(regions),
                                   '-o',str(root/'cds.png')],cwd=root,env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with (root/'cds.gene_pairs.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['passed'],'True')
            self.assertGreater((root/'cds.png').stat().st_size,1000)


if __name__=='__main__': unittest.main()
