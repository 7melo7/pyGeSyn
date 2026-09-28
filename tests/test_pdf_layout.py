from pathlib import Path
import re
import tempfile
import unittest

from pyGeSyn.plot import draw_synteny


def example(count, name='Genome', conn_ratio=1.5, output=None):
    regions = [(f'{name}_{i:03d}', 'chr1', 1, 2000, '') for i in range(count)]
    features = [[dict(category='te', strand='+', start=200, end=400)]
                for _ in regions]
    hits = [[dict(qstart=200, qend=800, sstart=200, send=800),
             dict(qstart=1200, qend=1600, sstart=1600, send=1200)]
            for _ in range(count - 1)]
    draw_synteny(regions, features, hits, output, conn_ratio)


def page_size(path):
    match = re.search(rb'/MediaBox\s*\[\s*0\s+0\s+([\d.]+)\s+([\d.]+)\s*\]',
                      path.read_bytes())
    if match is None:
        raise AssertionError('PDF page size not found')
    return tuple(map(float, match.groups()))


class PdfLayoutTests(unittest.TestCase):
    def test_page_grows_with_tracks_labels_and_spacing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sizes = {}
            for name, count, label, gap in (
                ('small', 2, 'Genome', 1.5),
                ('many', 40, 'Genome', 1.5),
                ('wide', 2, 'Very_long_genome_accession_' * 6, 1.5),
                ('spaced', 2, 'Genome', 8),
                ('single', 1, 'Genome', 1.5),
            ):
                path = root / f'{name}.pdf'
                example(count, label, gap, path)
                sizes[name] = page_size(path)
            self.assertGreater(sizes['many'][1], 11.69 * 72)
            self.assertGreater(sizes['many'][1], 8 * sizes['small'][1])
            self.assertGreater(sizes['wide'][0], sizes['small'][0] + 100)
            self.assertGreater(sizes['spaced'][1], sizes['small'][1])
            self.assertGreater(sizes['single'][1], 0)

    def test_invalid_layout(self):
        with self.assertRaises(ValueError):
            draw_synteny([], [], [], 'unused.pdf')
        with self.assertRaises(ValueError):
            example(2, conn_ratio=-1, output='unused.pdf')


if __name__ == '__main__':
    unittest.main()
