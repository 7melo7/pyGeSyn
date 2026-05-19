import argparse
import sys
from pathlib import Path

from .config import load_config
from .genome import extract_sequence
from .gff import parse_gene_structures, parse_te_features
from .blast import run_pairwise_blast
from .plot import draw_synteny
from .find import discover_regions, write_results
from .dedup import dedup_regions, write_dedup_results
from .dotplot import draw_dotplot


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


def cmd_plot(args):
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

    if args.sort:
        reverse = args.sort == 'desc'
        raw_regions.sort(key=lambda r: r[3] - r[2], reverse=reverse)

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


def cmd_find(args):
    workdir = Path("temp")
    workdir.mkdir(exist_ok=True)
    print("Discovering collinear regions ...", file=sys.stderr)
    all_candidates, best_per_genome = discover_regions(
        args.query, args.config,
        min_length=args.min_length,
        min_identity=args.min_identity,
        cluster_gap=args.cluster_gap,
        merge_gap=args.merge_gap,
        min_coverage=args.min_coverage,
        min_hits=args.min_hits,
        window_mult=args.window_mult,
        workdir=workdir,
    )

    query_region = None
    genome_order = None
    if '$' in args.query and not args.query.startswith('>'):
        parts = args.query.split('$')
        if len(parts) == 4:
            query_region = (parts[0], parts[1], parts[2], parts[3])
            config = load_config(args.config)
            genome_order = list(config.keys())

    write_results(all_candidates, best_per_genome,
                  args.all_output, args.best_output,
                  query_region=query_region, genome_order=genome_order)


def cmd_dedup(args):
    print("De-duplicating regions ...", file=sys.stderr)
    hap_info, nonredundant = dedup_regions(
        args.regions, args.config, args.coverage)
    write_dedup_results(hap_info, nonredundant,
                        args.hap_output, args.nr_output)


def cmd_dotplot(args):
    print("Drawing pairwise dotplot ...", file=sys.stderr)
    draw_dotplot(args.regions, args.config, args.pair, args.output,
                 args.min_length, args.min_identity)


def main():
    parser = argparse.ArgumentParser(
        description='pyGeSyn - Synteny visualization for genomic regions',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command')

    p_plot = sub.add_parser('plot', help='Draw synteny plot from regions file',
                            formatter_class=argparse.RawDescriptionHelpFormatter,
                            epilog='''\
Example:
  pyGeSyn plot genomes.json regions.csv -o synteny.pdf

Config file format (JSON):
{
  "GenomeA": {
    "fasta": "/path/to/genome.fa",
    "gff3": "/path/to/genes.gff3",
    "te": "/path/to/te.gff3"
  }
}

Regions file format (CSV):
  GenomeName,Chromosome,Start,End''')
    p_plot.add_argument('config', help='Config JSON file')
    p_plot.add_argument('regions', help='Regions CSV file')
    p_plot.add_argument('-o', '--output', default='synteny.png',
                        help='Output image (default: synteny.png)')
    p_plot.add_argument('--min-length', type=int, default=200)
    p_plot.add_argument('--min-identity', type=float, default=80.0)
    p_plot.add_argument('--conn-ratio', type=float, default=1.5)
    p_plot.add_argument('--sort', choices=['asc', 'desc'], default=None,
                        help='Sort regions by size: asc (smallest first) or desc (largest first)')
    p_plot.set_defaults(func=cmd_plot)

    p_find = sub.add_parser('find', help='Find collinear regions from query',
                            formatter_class=argparse.RawDescriptionHelpFormatter,
                            epilog='''\
Examples:
  pyGeSyn find GenomeA$Chr1$100000$150000 genomes.json
  pyGeSyn find query.fasta genomes.json
  pyGeSyn find ">query\\nATCG..." genomes.json''')
    p_find.add_argument('query', help='FASTA file, FASTA string, or Genome$Chr$Start$End')
    p_find.add_argument('config', help='Config JSON file')
    p_find.add_argument('--all-output', default='all_candidates.tsv')
    p_find.add_argument('--best-output', default='best_regions.csv')
    p_find.add_argument('--min-length', type=int, default=200)
    p_find.add_argument('--min-identity', type=float, default=80.0)
    p_find.add_argument('--cluster-gap', type=int, default=50000)
    p_find.add_argument('--merge-gap', type=int, default=200000)
    p_find.add_argument('--min-coverage', type=float, default=0.0)
    p_find.add_argument('--min-hits', type=int, default=2)
    p_find.add_argument('--window-mult', type=float, default=5,
                        help='Sliding window size as multiple of query length (default: 5)')
    p_find.set_defaults(func=cmd_find)

    p_dedup = sub.add_parser('dedup', help='Remove redundant regions by mutual coverage',
                             formatter_class=argparse.RawDescriptionHelpFormatter)
    p_dedup.add_argument('regions', help='Regions CSV file')
    p_dedup.add_argument('config', help='Config JSON file')
    p_dedup.add_argument('--coverage', type=float, default=0.8,
                         help='Mutual coverage threshold (default: 0.8)')
    p_dedup.add_argument('--hap-output', default='haplotypes.tsv')
    p_dedup.add_argument('--nr-output', default='nonredundant.csv')
    p_dedup.set_defaults(func=cmd_dedup)

    p_dot = sub.add_parser('dotplot', help='Pairwise diagonal collinearity plot',
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    p_dot.add_argument('regions', help='Regions CSV file')
    p_dot.add_argument('config', help='Config JSON file')
    p_dot.add_argument('pair', help='Two genome names: genomeA,genomeB')
    p_dot.add_argument('-o', '--output', default='dotplot.png')
    p_dot.add_argument('--min-length', type=int, default=100)
    p_dot.add_argument('--min-identity', type=float, default=80.0)
    p_dot.set_defaults(func=cmd_dotplot)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    args.func(args)
