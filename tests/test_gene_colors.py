import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from matplotlib.colors import to_rgba
from matplotlib.figure import Figure

from pyGeSyn.cli import main
from pyGeSyn.config import load_config, load_gene_colors
from pyGeSyn.gff import parse_gene_structures
from pyGeSyn.plot import _draw_gene, COLOR_GENE


GFF = '''chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene1;
chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=gene1.t1;Parent=gene1;
chr1\t.\texon\t1\t100\t.\t+\t.\tID=ex1;Parent=gene1.t1;
chr1\t.\tCDS\t21\t80\t.\t+\t0\tID=cds1;Parent=gene1.t1;
chr1\t.\texon\t200\t300\t.\t+\t.\tID=ex2;Parent=gene1.t1;
chr1\t.\tmRNA\t600\t1000\t.\t+\t.\tID=gene1.t2;Parent=gene1;
chr1\t.\tgene\t1100\t1400\t.\t-\t.\tID=gene10;
chr1\t.\tmRNA\t1100\t1400\t.\t-\t.\tID=gene10.t1;Parent=gene10;
chr1\t.\texon\t1100\t1400\t.\t-\t.\tID=ex3;Parent=gene10.t1;
'''


class GeneColorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.gff = self.root / 'genes.gff3'
        self.gff.write_text(GFF)
        self.colors = self.root / 'colors.txt'

    def genes(self, colors=None, end=1500):
        return parse_gene_structures(self.gff, 'chr1', 1, end, colors)

    def test_all_delimiters_and_hex_validation(self):
        self.colors.write_text('\ufeffgene1,#ff0000\ngene1.t2 #00ff00\ngene10\t#0000FF\n')
        self.assertEqual(load_gene_colors(self.colors),
                         {'gene1': '#FF0000', 'gene1.t2': '#00FF00', 'gene10': '#0000FF'})
        for text in ('gene1 red', 'gene1 #FFF', 'gene1 #GG0000', 'gene1 #123456 extra',
                     'gene1 #FF0000\ngene1 #00FF00'):
            self.colors.write_text(text)
            with self.assertRaises(ValueError):
                load_gene_colors(self.colors)

    def test_exact_gene_and_nonselected_transcript(self):
        genes = self.genes({'gene1': '#FF0000'})
        self.assertEqual(genes[0]['color'], '#FF0000')
        self.assertNotIn('color', genes[1])
        self.assertEqual(self.genes({'gene1.t2': '#00FF00'})[0]['color'], '#00FF00')
        # Transcript 2 lies outside the plotted region but belongs to gene1.
        self.assertEqual(self.genes({'gene1.t2': '#00FF00'}, end=400)[0]['color'], '#00FF00')
        self.assertNotIn('color', self.genes({'gene': '#FF0000'})[0])

    def test_gene_precedence_and_transcript_conflicts(self):
        colors = {'gene1.t1': '#FF0000', 'gene1.t2': '#00FF00'}
        with self.assertRaisesRegex(ValueError, 'conflicting transcript colors'):
            self.genes(colors)
        colors['gene1'] = '#0000FF'
        self.assertEqual(self.genes(colors)[0]['color'], '#0000FF')

    def test_exon_cds_utr_colors_and_default(self):
        genes = self.genes({'gene1.t2': '#FF0000'})
        for gene, expected in zip(genes, ('#FF0000', COLOR_GENE)):
            ax = Figure().add_subplot()
            _draw_gene(ax, gene, 0, 1, 20, 10)
            self.assertTrue(ax.patches)
            for rectangle in ax.patches:
                self.assertEqual(rectangle.get_facecolor(), to_rgba(expected))
            if ax.lines:
                self.assertEqual(ax.lines[0].get_color(), 'k')

    def test_optional_config_validation_and_duplicate_genomes(self):
        fasta = self.root / 'genome.fa'
        fasta.write_text('>chr1\nACGT')
        cfg = self.root / 'genomes.json'
        paths = {'fasta': str(fasta), 'gff3': str(self.gff)}
        cfg.write_text(json.dumps({'A': paths}))
        self.assertNotIn('gene_colors', load_config(str(cfg))['A'])
        paths['gene_colors'] = str(self.colors)
        cfg.write_text(json.dumps({'A': paths}))
        with self.assertRaises(FileNotFoundError):
            load_config(str(cfg))
        self.colors.write_text('gene1 #FF0000')
        self.assertEqual(load_config(str(cfg))['A']['gene_colors'], str(self.colors))
        cfg.write_text(json.dumps({'A': {'fasta': str(fasta),
                                        'gene_colors': str(self.colors)}}))
        with self.assertRaisesRegex(ValueError, 'gene_colors requires gff3'):
            load_config(str(cfg))
        cfg.write_text(json.dumps({'A': paths}))
        cfg2 = self.root / 'other.json'
        cfg2.write_text(json.dumps({'A': {'fasta': str(fasta), 'gff3': str(self.gff)}}))
        with self.assertRaisesRegex(ValueError, 'Duplicate genome'):
            load_config(f'{cfg},{cfg2}')

    def test_plot_cli_reads_config_and_applies_transcript_color(self):
        fasta = self.root / 'genome.fa'
        fasta.write_text('>chr1\n' + 'ACGT' * 400)
        self.colors.write_text('gene1.t2,#FF0000')
        cfg = self.root / 'genomes.json'
        cfg.write_text(json.dumps({'A': {'fasta': str(fasta), 'gff3': str(self.gff),
                                        'gene_colors': str(self.colors)}}))
        regions = self.root / 'regions.csv'
        regions.write_text('A,chr1,1,1500\n')
        with patch('sys.argv', ['pyGeSyn', 'plot', str(cfg), str(regions)]), \
             patch('pyGeSyn.cli._init_workdir', return_value=self.root), \
             patch('pyGeSyn.cli.draw_synteny') as draw:
            main()
        features = draw.call_args.args[1][0]
        self.assertEqual(features[0]['color'], '#FF0000')
        self.assertNotIn('color', features[1])


if __name__ == '__main__':
    unittest.main()
