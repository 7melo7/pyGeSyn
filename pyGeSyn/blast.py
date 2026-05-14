import subprocess
import os
import sys
from pathlib import Path


def run_blastn(seq_a, name_a, seq_b, name_b, min_length=200, min_identity=80.0,
               workdir=None):
    if workdir:
        fa_a = Path(workdir) / f"{name_a}.fa"
        fa_b = Path(workdir) / f"{name_b}.fa"
        out = Path(workdir) / f"{name_a}_vs_{name_b}_blast.txt"
        fa_a.write_text(f'>{name_a}\n{seq_a}\n', encoding='utf-8')
        fa_b.write_text(f'>{name_b}\n{seq_b}\n', encoding='utf-8')
        result = subprocess.run(
            ['blastn', '-query', str(fa_a), '-subject', str(fa_b),
             '-outfmt', '6', '-out', str(out)],
            capture_output=True, text=True
        )
        blast_output = out.read_text()
    else:
        import tempfile
        tmp_a = tmp_b = None
        try:
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.fa', delete=False, encoding='utf-8'
            ) as f:
                f.write(f'>{name_a}\n{seq_a}\n')
                tmp_a = f.name
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.fa', delete=False, encoding='utf-8'
            ) as f:
                f.write(f'>{name_b}\n{seq_b}\n')
                tmp_b = f.name
            result = subprocess.run(
                ['blastn', '-query', tmp_a, '-subject', tmp_b, '-outfmt', '6'],
                capture_output=True, text=True
            )
            blast_output = result.stdout
        finally:
            if tmp_a and os.path.exists(tmp_a):
                os.unlink(tmp_a)
            if tmp_b and os.path.exists(tmp_b):
                os.unlink(tmp_b)

    if result.returncode not in (0, 1):
        raise RuntimeError(f"blastn failed: {result.stderr}")

    hits = []
    for line in blast_output.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 12:
            continue

        pident = float(parts[2])
        length = int(parts[3])
        qstart = int(parts[6])
        qend = int(parts[7])
        sstart = int(parts[8])
        send = int(parts[9])

        if length < min_length or pident < min_identity:
            continue

        hits.append({
            'qstart': qstart,
            'qend': qend,
            'sstart': sstart,
            'send': send,
            'identity': pident,
            'length': length,
        })

    return hits


def run_pairwise_blast(regions, min_length=200, min_identity=80.0,
                       workdir=None):
    results = []
    n = len(regions)

    for i in range(n - 1):
        name_a, chrom_a, start_a, end_a, seq_a = regions[i]
        name_b, chrom_b, start_b, end_b, seq_b = regions[i + 1]

        label_a = f"{name_a}_{chrom_a}:{start_a}-{end_a}"
        label_b = f"{name_b}_{chrom_b}:{start_b}-{end_b}"

        safe_a = f"{name_a}_{chrom_a}_{start_a}_{end_a}"
        safe_b = f"{name_b}_{chrom_b}_{start_b}_{end_b}"

        print(f"  blastn: {name_a} vs {name_b} ...", file=sys.stderr, end=' ')
        try:
            hits = run_blastn(
                seq_a, safe_a, seq_b, safe_b,
                min_length, min_identity, workdir=workdir
            )
        except FileNotFoundError:
            raise RuntimeError(
                "blastn not found. Please ensure BLAST+ is installed and in PATH."
            )
        print(f"{len(hits)} hits", file=sys.stderr)
        results.append(hits)

    return results
