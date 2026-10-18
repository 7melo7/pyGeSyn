import subprocess
import sys
import os
import math
import tempfile
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left
from concurrent.futures import ThreadPoolExecutor, as_completed

from .genome import extract_sequence
from .config import load_config

MIN_LENGTH = 200
MIN_IDENTITY = 80.0


def discover_regions(query_input, config_path, min_length=MIN_LENGTH,
                     min_identity=MIN_IDENTITY, cluster_gap=None, merge_gap=None,
                     min_coverage=0.0, min_hits=2, window_mult=5,
                     workdir=None, threads=1, repeatmaskmode=0):
    if cluster_gap is not None or merge_gap is not None:
        raise ValueError('cluster_gap and merge_gap are unsupported by the sliding-window '
                         'algorithm; use window_mult instead')
    _validate_search(min_hits, window_mult, min_coverage)
    if threads < 1:
        raise ValueError('threads must be at least 1')
    if repeatmaskmode not in (0, 1):
        raise ValueError("repeatmaskmode must be 0 or 1")
    config = load_config(config_path)

    query_name, query_seq = _resolve_query(query_input, config)
    query_len = len(query_seq)
    if query_len == 0:
        raise ValueError('Query sequence must not be empty')

    print(f"Query: {query_name}  length={query_len}bp", file=sys.stderr)

    query_te, query_offset = None, 0
    if repeatmaskmode:
        query_te, query_offset = _query_te_context(query_input, config)

    genomes = [(g, p) for g, p in config.items()
               if not query_input.startswith(g + '$')]

    if threads <= 1 or len(genomes) <= 1:
        return _run_sequential(query_seq, query_len, genomes, min_length,
                               min_identity, min_coverage, min_hits,
                               window_mult, workdir, repeatmaskmode, query_te, query_offset)

    all_results = []
    best_per_genome = []

    def process_genome(genome_name, paths):
        blast_out = None
        if workdir:
            blast_out = Path(workdir) / f"blast_{genome_name}.txt"
        raw = _blast_genome(query_seq, paths['fasta'], min_length, min_identity,
                            blast_out)
        if not raw:
            return genome_name, [], None
        filtered = _apply_repeatmask(raw, paths, repeatmaskmode, query_te, query_offset)
        if not filtered:
            return genome_name, [], None
        candidates = _find_candidates(filtered, query_len, min_coverage, window_mult, min_hits)
        best = _pick_best(candidates) if candidates else None
        for c in candidates:
            c['genome'] = genome_name
        if best:
            best['genome'] = genome_name
        return genome_name, candidates, best

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_genome, g, p): g for g, p in genomes}
        for future in as_completed(futures):
            genome_name, candidates, best = future.result()
            print(f"  {genome_name}: {len(candidates)} candidates", file=sys.stderr)
            all_results.extend(candidates)
            if best:
                best_per_genome.append(best)

    return all_results, best_per_genome


def _run_sequential(query_seq, query_len, genomes, min_length, min_identity,
                    min_coverage, min_hits, window_mult, workdir,
                    repeatmaskmode=0, query_te=None, query_offset=0):
    all_results = []
    best_per_genome = []

    for genome_name, paths in genomes:
        print(f"  {genome_name}: blastn ...", file=sys.stderr, end=' ')
        blast_out = None
        if workdir:
            blast_out = Path(workdir) / f"blast_{genome_name}.txt"
        raw = _blast_genome(query_seq, paths['fasta'], min_length, min_identity,
                            blast_out)
        print(f"{len(raw)} raw", file=sys.stderr, end='')

        if not raw:
            print()
            continue

        filtered = _apply_repeatmask(raw, paths, repeatmaskmode, query_te, query_offset)
        print(f" -> {len(filtered)} filtered", file=sys.stderr, end='')

        if not filtered:
            print()
            continue

        candidates = _find_candidates(filtered, query_len, min_coverage, window_mult, min_hits)
        print(f" -> {len(candidates)} candidates", file=sys.stderr)

        for c in candidates:
            c['genome'] = genome_name
            all_results.append(c)

        best = _pick_best(candidates)
        if best:
            best['genome'] = genome_name
            best_per_genome.append(best)

    return all_results, best_per_genome


def _resolve_query(query_input, config):
    if os.path.isfile(query_input):
        seq = _parse_fasta(Path(query_input).read_text())
        return 'query', seq
    if query_input.startswith('>') or '\n' in query_input:
        return 'query', _parse_fasta(query_input)
    parts = query_input.split('$')
    if len(parts) == 4:
        genome, chrom, s, e = parts
        if genome not in config:
            raise ValueError(f"Genome '{genome}' not found")
        seq = extract_sequence(config[genome]['fasta'], chrom, int(s), int(e))
        return f"{genome}_{chrom}:{s}-{e}", seq
    raise ValueError("Query: FASTA file/string or Genome$Chr$Start$End")


def _parse_fasta(text):
    return ''.join(l.strip() for l in text.split('\n') if not l.startswith('>'))


def _blast_genome(query_seq, subject_fasta, min_length, min_identity,
                  blast_out=None):
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, encoding='utf-8'
    ) as f:
        f.write(f'>query\n{query_seq}\n')
        qf = f.name
    try:
        r = subprocess.run(
            ['blastn', '-query', qf, '-subject', subject_fasta,
             '-outfmt', '6 qseqid sseqid pident length '
                        'qstart qend sstart send evalue bitscore',
             '-evalue', '1e-5'],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"blastn: {r.stderr}")
        if blast_out:
            blast_out.parent.mkdir(parents=True, exist_ok=True)
            blast_out.write_text(r.stdout)
        hsps = []
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            p = line.split('\t')
            if len(p) < 10:
                continue
            pid = float(p[2])
            ln = int(p[3])
            qs, qe = int(p[4]), int(p[5])
            ss, se = int(p[6]), int(p[7])
            bs = float(p[9])
            if ln < min_length or pid < min_identity:
                continue
            if qs > qe:
                qs, qe = qe, qs
                ss, se = se, ss
            strand = '+' if ss <= se else '-'
            if ss > se:
                ss, se = se, ss
            hsps.append({
                'chrom': p[1], 'qstart': qs, 'qend': qe,
                'sstart': ss, 'send': se, 'strand': strand,
                'identity': pid, 'length': ln, 'bitscore': bs,
            })
        return hsps
    finally:
        os.unlink(qf)


def _load_te_index(path):
    """Index individual GFF intervals without merging neighboring TEs."""
    intervals = defaultdict(list)
    with open(path) as stream:
        for lineno, line in enumerate(stream, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{lineno}: expected nine GFF columns')
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f'{path}:{lineno}: invalid TE coordinates') from exc
            if start < 1 or end < start:
                raise ValueError(f'{path}:{lineno}: invalid TE interval')
            intervals[fields[0]].append((start, end))
    index = {}
    for chrom, regions in intervals.items():
        starts, max_ends = [], []
        furthest = 0
        for start, end in sorted(regions):
            starts.append(start)
            furthest = max(furthest, end)
            max_ends.append(furthest)
        index[chrom] = (starts, max_ends)
    return index


def _strictly_in_te(intervals, start, end):
    if not intervals:
        return False
    low, high = sorted((start, end))
    starts, max_ends = intervals
    # Only TEs starting strictly before low can contain this HSP.
    pos = bisect_left(starts, low) - 1
    return pos >= 0 and max_ends[pos] > high


def _query_te_context(query_input, config):
    if not os.path.isfile(query_input) and not query_input.startswith('>') and '\n' not in query_input:
        parts = query_input.split('$')
        if len(parts) == 4:
            genome, chrom, start, _ = parts
            path = config[genome].get('te')
            if path:
                return _load_te_index(path).get(chrom), max(1, int(start)) - 1
            print(f'  Query {genome}: no TE annotation; query-side masking skipped', file=sys.stderr)
            return None, 0
    print('  FASTA query has no genome coordinates; query-side masking skipped', file=sys.stderr)
    return None, 0


def _apply_repeatmask(hsps, paths, mode, query_te, query_offset):
    if not mode:
        return hsps
    path = paths.get('te')
    subject_te = _load_te_index(path) if path else {}
    if not path:
        print(f"  Subject {paths['fasta']}: no TE annotation; subject-side masking skipped",
              file=sys.stderr)
    return [h for h in hsps
            if not (_strictly_in_te(query_te, h['qstart'] + query_offset,
                                    h['qend'] + query_offset)
                    or _strictly_in_te(subject_te.get(h['chrom']), h['sstart'], h['send']))]


def _validate_search(min_hits, window_mult, min_coverage):
    if not isinstance(min_hits, int) or min_hits < 1:
        raise ValueError('min_hits must be a positive integer')
    if not math.isfinite(window_mult) or window_mult <= 0:
        raise ValueError('window_mult must be finite and greater than zero')
    if not 0 <= min_coverage <= 1:
        raise ValueError('min_coverage must be between 0 and 1')


def _find_candidates(hsps, query_len, min_coverage, window_mult=5, min_hits=2):
    _validate_search(min_hits, window_mult, min_coverage)
    if query_len <= 0:
        raise ValueError('Query sequence must not be empty')
    by_cs = defaultdict(list)
    for h in hsps:
        by_cs[h['chrom']].append(h)
    candidates = []
    for chrom, group in by_cs.items():
        if len(group) < min_hits:
            continue
        group.sort(key=lambda h: h['sstart'])
        window = max(1, int(query_len * window_mult))
        step = max(1, window // 2)
        buffer = window // 4
        s_min_all = group[0]['sstart']
        s_max_all = max(h['send'] for h in group)
        w_start = s_min_all
        while w_start <= s_max_all:
            w_end = w_start + window
            seed = [h for h in group if w_start <= h['sstart'] <= w_end]
            if len(seed) >= min_hits:
                q_low = min(h['qstart'] for h in seed)
                q_high = max(h['qend'] for h in seed)
                s_low = min(h['sstart'] for h in seed)
                s_high = max(h['send'] for h in seed)
                region = [h for h in group
                          if s_low - buffer <= h['sstart'] <= s_high + buffer]
                q_union = _union(region, 'qstart', 'qend')
                coverage = q_union / query_len
                if coverage >= min_coverage:
                    ordered = sorted(region, key=lambda h: h['qstart'])
                    coll = _collinearity(ordered)
                    score = coverage * (0.8 + 0.2 * coll)

                    near_start = [h for h in ordered
                                  if h['qstart'] <= query_len * 0.10]
                    near_end = [h for h in ordered
                                if h['qend'] >= query_len * 0.90]

                    if near_start:
                        best_start = max(near_start, key=lambda h: h['length'])
                    else:
                        best_start = ordered[0]
                    if near_end:
                        best_end = max(near_end, key=lambda h: h['length'])
                    else:
                        best_end = ordered[-1]

                    # Enclose all contributing HSPs, including reverse or mixed hits.
                    s_start = min(h['sstart'] for h in region)
                    s_end = max(h['send'] for h in region)
                    strands = {h['strand'] for h in region}
                    orientation = next(iter(strands)) if len(strands) == 1 else 'mixed'
                    total_length = sum(h['length'] for h in region)
                    identity = sum(h['identity'] * h['length'] for h in region) / total_length
                    q_nearest_start = best_start['qstart']
                    q_nearest_end = best_end['qend']
                    len_start = best_start['length']
                    len_end = best_end['length']
                    q_start_pct = q_nearest_start / query_len
                    q_end_pct = 1.0 - q_nearest_end / query_len
                    edge_bonus = max(0, 1.0 - q_start_pct - q_end_pct) * 0.15
                    len_bonus = min(1.0, (len_start + len_end) / 2000.0) * 0.10
                    score += edge_bonus + len_bonus

                    candidates.append({
                        'chrom': chrom,
                        'start': s_start,
                        'end': s_end,
                        'coverage': coverage,
                        'identity': identity,
                        'score': score,
                        'orientation': orientation,
                        'num_hits': len(ordered),
                        'subject_span': s_end - s_start + 1,
                        '_q_low': min(h['qstart'] for h in region),
                        '_q_high': max(h['qend'] for h in region),
                    })
            w_start += step
    candidates.sort(key=lambda c: c['score'], reverse=True)
    return candidates


def _collinearity(ordered):
    if len(ordered) < 2:
        return 1.0
    diffs = []
    for k in range(len(ordered) - 1):
        h1, h2 = ordered[k], ordered[k + 1]
        dq = h2['qstart'] - h1['qstart']
        s1 = -h1['sstart'] if h1['strand'] == '-' else h1['sstart']
        s2 = -h2['sstart'] if h2['strand'] == '-' else h2['sstart']
        ds = s2 - s1
        if dq > 0:
            diffs.append(abs(dq - ds) / max(dq, abs(ds)))
        else:
            diffs.append(1.0)
    avg = sum(diffs) / len(diffs)
    return 1.0 / (1.0 + avg)


def _union(items, sk, ek):
    if not items:
        return 0
    iv = sorted((it[sk], it[ek]) for it in items)
    m = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return sum(e - s + 1 for s, e in m)


def _pick_best(candidates):
    return max(candidates, key=lambda c: c['score']) if candidates else None


def write_results(all_candidates, best_per_genome, all_path, best_path,
                   query_region=None, genome_order=None):
    with open(all_path, 'w') as f:
        f.write("genome\tchrom\tstart\tend\tcoverage\tidentity\tscore\t"
                "orientation\tnum_hits\tsubject_span\tq_low\tq_high\n")
        for c in all_candidates:
            f.write(f"{c['genome']}\t{c['chrom']}\t{c['start']}\t{c['end']}\t"
                    f"{c['coverage']:.4f}\t{c['identity']:.2f}\t{c['score']:.4f}\t"
                    f"{c['orientation']}\t{c['num_hits']}\t{c['subject_span']}\t"
                    f"{c.get('_q_low', 0)}\t{c.get('_q_high', 0)}\n")
    with open(best_path, 'w') as f:
        f.write("# Best collinear region per genome\n")
        f.write("# (regions.csv format, ready for pyGeSyn plot)\n")
        if query_region:
            f.write(f"{query_region[0]},{query_region[1]},"
                    f"{query_region[2]},{query_region[3]}\n")
        best_map = {c['genome']: c for c in best_per_genome}
        order = genome_order or []
        for name in order:
            if name in best_map:
                c = best_map[name]
                f.write(f"{c['genome']},{c['chrom']},{c['start']},{c['end']}\n")
        for c in sorted(best_per_genome, key=lambda x: x['score'], reverse=True):
            if c['genome'] not in order:
                f.write(f"{c['genome']},{c['chrom']},{c['start']},{c['end']}\n")
    print(f"  All candidates -> {all_path}", file=sys.stderr)
    print(f"  Best per genome -> {best_path}", file=sys.stderr)
