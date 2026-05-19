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
pyGeSyn find MSU$Chr3$35647669$35741607 genomes.json

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
| `--min-length` | `100` | Minimum HSP length (bp) |
| `--min-identity` | `70` | Minimum HSP identity (%) |
| `--min-coverage` | `0.0` | Minimum query coverage (0–1) |
| `--min-hits` | `2` | Minimum HSPs per candidate |
| `--window-mult` | `5` | Sliding window as multiple of query length |
| `--threads` | `1` | Parallel threads for multi-genome search |

**Algorithm**: blastn → filter repeats (chromosome-aware) → group by chromosome → sliding window (5× query, 50% overlap) → expand & score by coverage + collinearity.

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
        "te": "/data/MSU.te.gff3"
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
- A4 output, auto-detect system font (Arial / Helvetica / DejaVu Sans)

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

`find` saves raw blast results to `temp/blast_*.txt`.
`plot` saves region sequences and pairwise blast to `temp/`.
`dedup` uses temporary files (auto-cleaned).
