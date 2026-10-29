import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from matplotlib.figure import Figure
from matplotlib.colors import to_rgba
from pyGeSyn import dotplot
from pyGeSyn.gff import parse_gene_structures
from test_gene_colors import GFF


class DotplotColorTests(unittest.TestCase):
    def test_both_tracks_custom_and_default_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp)/'genes.gff'
            gff.write_text(GFF)
            genes = parse_gene_structures(gff,'chr1',1,1500,{'gene1.t2':'#FF0000'})
            for draw in (dotplot._gene_h,dotplot._gene_v):
                for gene,color in zip(genes,('#FF0000',dotplot.COLOR_GENE)):
                    ax=Figure().add_subplot()
                    draw(ax,gene,0,1)
                    self.assertTrue(ax.patches)
                    self.assertTrue(all(p.get_facecolor()==to_rgba(color) for p in ax.patches))
                    self.assertTrue(all(line.get_color()=='k' for line in ax.lines))

    def test_config_colors_follow_genomes_after_axis_swap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            gff=root/'genes.gff';gff.write_text(GFF)
            fa=root/'genome.fa';fa.write_text('>chr1\n'+'A'*2000)
            ca=root/'a.txt';ca.write_text('gene1.t2 #FF0000\n')
            cb=root/'b.txt';cb.write_text('gene1,#00FF00\n')
            cfg=root/'genomes.json'
            cfg.write_text(json.dumps({name:dict(fasta=str(fa),gff3=str(gff),gene_colors=str(color))
                                       for name,color in [('A',ca),('B',cb)]}))
            regions=root/'regions.csv'
            for swapped in (False,True):
                regions.write_text('A,chr1,1,1500\nB,chr1,1,2000\n' if swapped else
                                   'A,chr1,1,2000\nB,chr1,1,1500\n')
                with patch.object(dotplot,'create_workdir',return_value=root), \
                     patch.object(dotplot,'_blast_pair',return_value=[]), \
                     patch.object(dotplot,'_draw_track_h',wraps=dotplot._draw_track_h) as horizontal, \
                     patch.object(dotplot,'_draw_track_v',wraps=dotplot._draw_track_v) as vertical:
                    output=root/f'colors-{swapped}.png'
                    dotplot.draw_dotplot(regions,str(cfg),'A,B',output)
                    self.assertGreater(output.stat().st_size,1000)
                    self.assertEqual(horizontal.call_args.args[1][0]['color'], '#00FF00' if swapped else '#FF0000')
                    self.assertEqual(vertical.call_args.args[1][0]['color'], '#FF0000' if swapped else '#00FF00')


if __name__=='__main__':
    unittest.main()
