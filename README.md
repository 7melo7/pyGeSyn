# pyGeSyn

Synteny visualization tool for genomic regions — draw collinearity relationships between multiple genomes with gene/TE annotations.

## Requirements

- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (make sure `blastn` is in your `PATH`)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) / Anaconda

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/pyGeSyn.git
cd pyGeSyn

# 2. Create and activate the conda environment
conda create -n pyGeSyn python=3.12 matplotlib -y
conda activate pyGeSyn

# 3. Install the package
pip install -e .
```

Verify installation:
```bash
python -m pyGeSyn --help
```

## Configuration Files

### 1. Genome config (`genomes.json`)

JSON file mapping each genome's unique name to its file paths:

```json
{
    "Osa": {
        "fasta": "/data/genomes/Osa.fasta",
        "gff3": "/data/annotations/Osa.gene.gff3",
        "te": "/data/annotations/Osa.te.gff3"
    },
    "Oru": {
        "fasta": "/data/genomes/Oru.fasta",
        "gff3": "/data/annotations/Oru.gene.gff3",
        "te": "/data/annotations/Oru.te.gff3"
    }
}
```

| Field  | Description                                              |
|--------|----------------------------------------------------------|
| fasta  | Genome FASTA sequence file                               |
| gff3   | Gene annotation GFF3 (must contain gene→mRNA→exon/CDS)   |
| te     | Transposable element annotation GFF3                     |

Paths can be absolute or relative to the working directory.

### 2. Regions file (`regions.csv`)

CSV file specifying the genomic regions to compare. One region per line, comma-separated:

```
GenomeName,Chromosome,Start,End
```

Example:
```
Osa,Chr1,40000,50000
Oru,Chr1,43000,60000
Opu,Chr1,52000,60000
```

Lines starting with `#` are treated as comments. Coordinates are 1-based inclusive (GFF3 convention).

Each adjacent pair of regions is compared via blastn (line1 vs line2, line2 vs line3, ...).

## Usage

```bash
conda activate pyGeSyn
pyGeSyn genomes.json regions.csv [OPTIONS]
```

### Options

| Option              | Default        | Description                                   |
|---------------------|----------------|-----------------------------------------------|
| `-o`, `--output`    | `synteny.png`  | Output image file (supports png/pdf/svg)      |
| `--min-length`      | `200`          | Minimum blast hit length (bp)                 |
| `--min-identity`    | `80`           | Minimum blast hit identity (%)                |
| `--conn-ratio`      | `1.5`          | Connection-area height to feature height ratio |

### Example

```bash
pyGeSyn genomes.json regions.csv -o synteny.pdf --min-length 300 --min-identity 85
```

## Output

- **Synteny plot** — A4-sized image showing:
  - Each genome region as a horizontal gray bar
  - Genes: blue rectangles for CDS (wide) / UTR (narrow), black lines for introns
  - TEs: orange rectangles
  - Positive-strand features drawn above the bar, negative-strand below
  - Collinearity connections (bezier ribbons) between adjacent tracks:
    - Gray: collinear alignments
    - Red: inverted alignments
  - Legend below the plot

- **Intermediate files** — saved in `temp/`:
  - `{name}_{chr}_{start}_{end}.fa` — extracted region sequences
  - `{nameA}_vs_{nameB}_blast.txt` — blastn output (outfmt 6)

## GFF3 Format Requirements

### Gene GFF3
Must contain the parent-child hierarchy:
```
Chr1  .  gene  1000  5000  .  +  .  ID=gene1
Chr1  .  mRNA  1000  5000  .  +  .  ID=mRNA1;Parent=gene1
Chr1  .  exon  1000  1200  .  +  .  Parent=mRNA1
Chr1  .  CDS   1050  1200  .  +  .  Parent=mRNA1
Chr1  .  exon  2000  2500  .  +  .  Parent=mRNA1
Chr1  .  CDS   2000  2450  .  +  .  Parent=mRNA1
```

- Only the **first mRNA** per gene is used
- UTR regions are automatically computed as exon regions not covered by CDS
- If no CDS annotations exist, the full exon is drawn

### TE GFF3
Any GFF3 features in the TE file are drawn as transposons. All feature types are accepted.
