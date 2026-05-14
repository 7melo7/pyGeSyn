import argparse
import sys
from pathlib import Path

from .config import load_config
from .genome import extract_sequence
from .gff import parse_gene_structures, parse_te_features
from .blast import run_pairwise_blast
from .plot import draw_synteny


def parse_regions_file(path):
    regions = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) != 4:
                raise ValueError(
                    f"Line {lineno}: expected 4 comma-separated fields, "
                    f"got {len(parts)}"
                )
            name, chrom, start_str, end_str = parts
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                raise ValueError(
                    f"Line {lineno}: start/end must be integers"
                )
            if start > end:
                raise ValueError(
                    f"Line {lineno}: start ({start}) > end ({end})"
                )
            regions.append((name, chrom, start, end))
    return regions


def main():
    parser = argparse.ArgumentParser(
        description='pyGeSyn - Synteny visualization for genomic regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''\
Example:
  pyGeSyn genomes.json regions.csv -o synteny.pdf --min-length 300 --min-identity 85

Config file format (JSON):
{
  "GenomeA": {
    "fasta": "/path/to/genome.fa",
    "gff3": "/path/to/genes.gff3",
    "te": "/path/to/te.gff3"
  },
  "GenomeB": { ... }
}

Regions file format (CSV, one line per region):
  GenomeName,Chromosome,Start,End
  Osa,Chr1,40000,50000
  Oru,Chr1,43000,60000
'''
    )
    parser.add_argument(
        'config', help='Configuration JSON file with genome paths'
    )
    parser.add_argument(
        'regions', help='Genome regions file (CSV: name,chrom,start,end)'
    )
    parser.add_argument(
        '-o', '--output', default='synteny.png',
        help='Output image file (default: synteny.png)'
    )
    parser.add_argument(
        '--min-length', type=int, default=200,
        help='Minimum blast hit length in bp (default: 200)'
    )
    parser.add_argument(
        '--min-identity', type=float, default=80.0,
        help='Minimum blast hit identity in percent (default: 80)'
    )
    parser.add_argument(
        '--conn-ratio', type=float, default=3,
        help='Ratio of connection-area height to feature height (default: 3)'
    )

    args = parser.parse_args()

    workdir = Path("temp")
    workdir.mkdir(exist_ok=True)

    print("[1/5] Loading configuration ...", file=sys.stderr)
    config = load_config(args.config)

    print("[2/5] Parsing regions file ...", file=sys.stderr)
    raw_regions = parse_regions_file(args.regions)

    if len(raw_regions) < 1:
        print("Error: at least one genome region is required", file=sys.stderr)
        sys.exit(1)

    for name, chrom, start, end in raw_regions:
        if name not in config:
            print(f"Error: genome '{name}' not found in config", file=sys.stderr)
            sys.exit(1)

    print("[3/5] Extracting sequences and annotations ...", file=sys.stderr)
    regions = []
    all_features = []

    for name, chrom, start, end in raw_regions:
        genome_cfg = config[name]
        seq = extract_sequence(genome_cfg['fasta'], chrom, start, end)

        fa_path = workdir / f"{name}_{chrom}_{start}_{end}.fa"
        fa_path.write_text(
            f">{name}_{chrom}:{start}-{end}\n{seq}\n", encoding='utf-8')

        genes = parse_gene_structures(genome_cfg['gff3'], chrom, start, end)
        tes = parse_te_features(genome_cfg['te'], chrom, start, end)

        combined = genes + tes
        combined.sort(key=lambda x: x.get('gene_start', x.get('start', 0)))

        regions.append((name, chrom, start, end, seq))
        all_features.append(combined)

        print(
            f"  {name} {chrom}:{start}-{end}  "
            f"seq={len(seq)}bp  genes={len(genes)}  TEs={len(tes)}",
            file=sys.stderr
        )

    print("[4/5] Running pairwise blastn ...", file=sys.stderr)
    blast_results = []
    if len(regions) >= 2:
        blast_results = run_pairwise_blast(
            regions, args.min_length, args.min_identity, workdir=workdir
        )

    print("[5/5] Generating synteny plot ...", file=sys.stderr)
    draw_synteny(regions, all_features, blast_results, args.output,
                 conn_ratio=args.conn_ratio)

    print(f"Done. Output: {args.output}", file=sys.stderr)
    print(f"  Intermediate files saved in: {workdir.resolve()}", file=sys.stderr)
