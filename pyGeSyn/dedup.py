import subprocess
import sys
import os
import tempfile
from pathlib import Path
from collections import defaultdict

from .genome import extract_sequence
from .config import load_config


def dedup_regions(regions_path, config_path, coverage_threshold=0.8):
    config = load_config(config_path)

    raw_regions = _parse_regions(regions_path)
    if len(raw_regions) < 2:
        print("Need at least 2 regions", file=sys.stderr)
        return [], []

    seqs = []
    labels = []
    for name, chrom, start, end in raw_regions:
        if name not in config:
            raise ValueError(f"Genome '{name}' not in config")
        seq = extract_sequence(config[name]['fasta'], chrom, start, end)
        seqs.append(seq)
        labels.append(f"{name}_{chrom}:{start}-{end}")
        print(f"  {name} {chrom}:{start}-{end}  {len(seq)}bp", file=sys.stderr)

    n = len(seqs)

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_path = os.path.join(tmpdir, 'all.fa')
        with open(fa_path, 'w') as f:
            for i, (lab, seq) in enumerate(zip(labels, seqs)):
                f.write(f'>{i}\n{seq}\n')

        print("  makeblastdb ...", file=sys.stderr)
        subprocess.run(
            ['makeblastdb', '-in', fa_path, '-dbtype', 'nucl', '-out',
             os.path.join(tmpdir, 'db')],
            capture_output=True, check=True)

        print("  blastn all-vs-all ...", file=sys.stderr)
        r = subprocess.run(
            ['blastn', '-query', fa_path, '-db', os.path.join(tmpdir, 'db'),
             '-outfmt', '6 qseqid sseqid pident length qstart qend'],
            capture_output=True, text=True)

        hits_by_pair = defaultdict(list)
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            p = line.split('\t')
            if len(p) < 7:
                continue
            qi = int(p[0])
            si = int(p[1])
            ln = int(p[3])
            qs, qe = int(p[4]), int(p[5])
            if qi == si:
                continue
            if qs > qe:
                qs, qe = qe, qs
            hits_by_pair[(qi, si)].append((qs, qe, ln))

        coverage = {}
        for (qi, si), hsps in hits_by_pair.items():
            cov = _union_len(hsps) / len(seqs[qi])
            coverage[(qi, si)] = cov

    mutual = {}
    for qi in range(n):
        for si in range(n):
            if qi == si:
                continue
            cov_ij = coverage.get((qi, si), 0.0)
            cov_ji = coverage.get((si, qi), 0.0)
            mutual[(qi, si)] = min(cov_ij, cov_ji)

    assigned = [False] * n
    haplotypes = []

    for i in range(n):
        if assigned[i]:
            continue
        group = [i]
        assigned[i] = True
        for j in range(n):
            if assigned[j]:
                continue
            if mutual.get((i, j), 0.0) >= coverage_threshold:
                group.append(j)
                assigned[j] = True
        haplotypes.append(group)

    hap_info = []
    nonredundant = []

    for idx, group in enumerate(haplotypes):
        names = [raw_regions[i][0] for i in group]
        hap_info.append((f"Hap{idx+1}", names))
        first = group[0]
        nonredundant.append(raw_regions[first])

    return hap_info, nonredundant


def _parse_regions(path):
    regions = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) != 4:
                continue
            regions.append((parts[0], parts[1], int(parts[2]), int(parts[3])))
    return regions


def _union_len(hsps):
    if not hsps:
        return 0
    intervals = sorted(hsps)
    merged = [list(intervals[0][:2])]
    for s, e, _ in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)


def write_dedup_results(hap_info, nonredundant, hap_path, nr_path):
    with open(hap_path, 'w') as f:
        for hap_id, names in hap_info:
            f.write(f"{hap_id}\t{', '.join(names)}\n")

    with open(nr_path, 'w') as f:
        f.write("# Non-redundant regions (regions.csv format)\n")
        for name, chrom, start, end in nonredundant:
            f.write(f"{name},{chrom},{start},{end}\n")

    print(f"  Haplotypes → {hap_path}", file=sys.stderr)
    print(f"  Non-redundant → {nr_path}", file=sys.stderr)
