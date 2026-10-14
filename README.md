# pyGeSyn

Synteny visualization toolkit — discover and draw collinearity relationships between genomic regions across multiple genomes.

## Requirements

- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (`blastn` in `PATH`)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Installation

```bash
git clone https://github.com/7melo7/pyGeSyn.git
cd pyGeSyn
conda create -n pyGeSyn python=3.12 matplotlib -y
conda activate pyGeSyn
pip install -e .
python -m pyGeSyn --help
```

## Quick Start

```bash
# 1. Find homologs: discover collinear regions in other genomes
pyGeSyn find 'MSU$Chr3$35647669$35741607' genomes.json

# 2. Plot: draw synteny from the generated regions file
pyGeSyn plot genomes.json best_regions.csv -o synteny.pdf
```

## Commands

### `pyGeSyn find` — Discover Homologous Regions

Finds the best collinear region in each target genome for a given query.

```
pyGeSyn find <query> <config> [options]
```

**query**: a FASTA file, inline FASTA string, or genome region `Genome$Chr$Start$End`

| Option | Default | Description |
|--------|---------|-------------|
| `--all-output` | `all_candidates.tsv` | All candidate regions with scores |
| `--best-output` | `best_regions.csv` | Best region per genome (ready for `plot`) |
| `--min-length` | `200` | Minimum HSP length (bp) |
| `--min-identity` | `80` | Minimum HSP identity (%) |
| `--min-coverage` | `0.0` | Minimum query coverage (0–1) |
| `--min-hits` | `2` | Minimum HSPs per candidate |
| `--window-mult` | `5` | Sliding window as multiple of query length |
| `--threads` | `1` | Parallel threads for multi-genome search |

**Algorithm**: blastn → filter repeats (chromosome-aware) → group by chromosome → sliding window (5× query, 50% overlap) → expand & score by coverage + collinearity.

Candidate coordinates enclose all contributing HSPs and always use ascending,
1-based inclusive coordinates. `subject_span` includes both endpoints.
`identity` is the alignment-length-weighted mean HSP identity (overlapping HSPs
are included independently); `orientation` is `+`, `-`, or `mixed`.
`--min-hits` controls the minimum number of HSPs in a seed window, including
single-HSP candidates when set to 1. `--window-mult` must be positive.
The previously ignored `--cluster-gap` and `--merge-gap` options have been
removed: the current algorithm uses sliding windows, not gap-based clustering.
Explicit Python API use of `cluster_gap` or `merge_gap` now raises an error.
Quote coordinate queries containing `$` to prevent shell variable expansion.

### `pyGeSyn plot` — Draw Synteny

Draws a collinearity plot from a regions CSV file.

```
pyGeSyn plot <config> <regions> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output` | `synteny.png` | Output image (png/pdf/svg) |
| `--min-length` | `200` | Minimum blast hit length (bp) |
| `--min-identity` | `80` | Minimum blast hit identity (%) |
| `--conn-ratio` | `1.5` | Connection area to feature height ratio |
| `--sort` | (none) | Sort regions by size: `asc` or `desc` |

### `pyGeSyn dedup` — Remove Redundant Regions

Groups similar/redundant genome regions into haplotypes by mutual coverage.

```
pyGeSyn dedup <regions> <config> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--coverage` | `0.8` | Mutual coverage threshold (0–1) |
| `--hap-output` | `haplotypes.tsv` | Haplotype group definitions |
| `--nr-output` | `nonredundant.csv` | Non-redundant regions for `plot` |

**Workflow**: extract sequences → `makeblastdb` → `blastn all-vs-all` → mutual coverage → greedy clustering. Each haplotype's first region becomes the representative.

### `pyGeSyn dotplot` — Pairwise Diagonal Synteny

Draws a 1-to-1 diagonal collinearity dotplot with gene/TE tracks on the axes.

```
pyGeSyn dotplot <regions> <config> <genomeA,genomeB> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output` | `dotplot.png` | Output image |
| `--min-length` | `100` | Minimum blast hit length (bp) |
| `--min-identity` | `80` | Minimum blast hit identity (%) |

The longer region is automatically placed on the X-axis.

## Configuration (`genomes.json`)

```json
{
    "MSU": {
        "fasta": "/data/MSU.fasta",
        "gff3": "/data/MSU.gene.gff3",
        "te": "/data/MSU.te.gff3",
        "gene_colors": "/data/MSU.gene_colors.txt"
    },
    "Oru": {
        "fasta": "/data/Oru.fasta",
        "gff3": "/data/Oru.gene.gff3",
        "te": "/data/Oru.te.gff3"
    }
}
```

| Field | Description |
|-------|-------------|
| `fasta` | Genome sequence (FASTA) |
| `gff3` | Gene annotation (must have gene→mRNA→exon/CDS hierarchy) |
| `te` | TE annotation (all feature types accepted) |
| `gene_colors` | Optional gene/transcript ID and color list, used by `plot` |

### Custom gene colors

Add `gene_colors` to any genome entry to color selected genes in `pyGeSyn plot`.
Omit it to keep the default gene color. Paths use the same working-directory
rules as `fasta`, `gff3`, and `te`. Gene coloring requires that genome's `gff3`.

The file has no header: each line contains an exact ID and a six-digit hexadecimal
RGB color (`#RRGGBB`), separated by a comma, space, or tab. For example:

```text
OsGene001,#FF0000
OsGene002 #00AA00
OsGene003.t2\t#6633CC
```

In the last example, replace `\t` with a real tab character. Blank lines and
lines starting with `#` are ignored. Hex digits are case-insensitive.

IDs match GFF3 attributes exactly: either `ID` on a `gene` record or `ID` on
an `mRNA`/`transcript` record whose `Parent` identifies that gene. Listing any
of a gene's transcripts colors the entire displayed gene, even when another
transcript is chosen for drawing or the listed transcript lies outside the
displayed subregion. IDs are not matched by prefix or by stripping suffixes.

All displayed exon, CDS, and UTR rectangles of a matching gene use its assigned
color; introns remain black, and unlisted genes keep the default blue. TE and
alignment colors are unchanged. Lists are independent for each genome.
Unknown IDs have no effect.

An explicit gene ID color takes precedence over transcript colors for that gene.
If only transcript IDs are listed, they must agree on the color; conflicting
colors raise an error. Invalid color codes and conflicting duplicate IDs also
raise an error with the file and line number.

```bash
pyGeSyn plot genomes.json regions.csv -o colored_synteny.pdf
```

## Regions Format (`regions.csv`)

```
GenomeName,Chromosome,Start,End
MSU,Chr3,35647669,35741607
Oru,Chr3,33771548,33822973
```

1-based inclusive coordinates. `#` for comments. Adjacent rows are blastn-compared pairwise.

## Plot Features

- **Genes**: CDS (wide blue) / UTR (narrow blue) / introns (black lines)
- **TEs**: orange rectangles
- **+ strand** above track bar, **− strand** below
- **Connections**: bezier ribbons, gray=collinear, red=inverted
- Automatic canvas height based on track count; export bounds include long
  labels and the legend with padding, so PDF pages are not limited to A4.
  PNG and SVG exports use the same content-fitting layout.
- Auto-detect system font (Arial / Helvetica / DejaVu Sans)

## Example

[synteny.pdf](synteny.pdf) — three rice genome regions (MSU, Oru, Oni) collinearity plot.

## GFF3 Requirements

Gene GFF3 must have parent-child hierarchy. Only the **first mRNA** per gene is used. UTR = exon − CDS (auto-computed).

```
Chr1  .  gene  1000  5000  .  +  .  ID=gene1
Chr1  .  mRNA  1000  5000  .  +  .  ID=mRNA1;Parent=gene1
Chr1  .  exon  1000  1200  .  +  .  Parent=mRNA1
Chr1  .  CDS   1050  1200  .  +  .  Parent=mRNA1
```

## Intermediate Files

`find`, `plot`, and `dotplot` each create a unique `temp/pyGeSyn-*` directory
and retain their intermediate files there. Existing files are never cleared.
`find` saves raw blast results as `blast_*.txt` inside that directory.
`plot` and `dotplot` save region sequences and pairwise blast results there.
Use distinct output filenames when running multiple commands concurrently.
`dedup` uses temporary files (auto-cleaned).

## Tests

Run `python -m unittest discover -s tests -v` in the installed environment.
Tests requiring BLAST+ are skipped when `blastn` is unavailable.
