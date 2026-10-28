import subprocess
import sys
import os
import math
import tempfile
from pathlib import Path
from collections import defaultdict
from bisect import bisect_right
from concurrent.futures import ThreadPoolExecutor, as_completed

from .genome import extract_sequence
from .config import load_config
from .chains import candidates as _chain_candidates

MIN_LENGTH = 200
MIN_IDENTITY = 80.0


def discover_regions(query_input, config_path, min_length=MIN_LENGTH,
                     min_identity=MIN_IDENTITY, cluster_gap=None, merge_gap=None,
                     min_coverage=0.0, min_hits=2, window_mult=5,
                     workdir=None, threads=1, repeatmaskmode=0, boundary_mode='chains'):
    if cluster_gap is not None or merge_gap is not None:
        raise ValueError('cluster_gap and merge_gap are unsupported by the sliding-window '
                         'algorithm; use window_mult instead')
    _validate_search(min_hits, window_mult, min_coverage)
    if threads < 1:
        raise ValueError('threads must be at least 1')
    _validate_repeatmaskmode(repeatmaskmode)
    if boundary_mode not in ('chains', 'anchors', 'window'):
        raise ValueError('boundary_mode must be chains, anchors, or window')
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
                               window_mult, workdir, repeatmaskmode, query_te, query_offset, boundary_mode)

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
        candidates = _find_candidates(filtered, query_len, min_coverage, window_mult, min_hits, boundary_mode)
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

    statuses = defaultdict(int)
    for candidate in best_per_genome:
        statuses[candidate.get('support_status', candidate.get('boundary_method', 'window'))] += 1
    print('  Selected support: ' + ', '.join(f'{key}={value}' for key, value in sorted(statuses.items())), file=sys.stderr)
    return all_results, best_per_genome


def _run_sequential(query_seq, query_len, genomes, min_length, min_identity,
                    min_coverage, min_hits, window_mult, workdir,
                    repeatmaskmode=0, query_te=None, query_offset=0, boundary_mode='chains'):
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

        candidates = _find_candidates(filtered, query_len, min_coverage, window_mult, min_hits, boundary_mode)
        print(f" -> {len(candidates)} candidates", file=sys.stderr)

        for c in candidates:
            c['genome'] = genome_name
            all_results.append(c)

        best = _pick_best(candidates)
        if best:
            best['genome'] = genome_name
            best_per_genome.append(best)

    statuses = defaultdict(int)
    for candidate in best_per_genome:
        statuses[candidate.get('support_status', candidate.get('boundary_method', 'window'))] += 1
    print('  Selected support: ' + ', '.join(f'{key}={value}' for key, value in sorted(statuses.items())), file=sys.stderr)
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
    """Index the TE union with prefix lengths for logarithmic overlap queries."""
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
        merged = []
        for start, end in sorted(regions):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        starts, ends, prefix = [], [], [0]
        for start, end in merged:
            starts.append(start)
            ends.append(end)
            prefix.append(prefix[-1] + end - start + 1)
        index[chrom] = (starts, ends, prefix)
    return index


def _te_coverage(intervals, start, end):
    """Fraction of a 1-based inclusive interval covered by the TE union."""
    if not intervals:
        return 0.0
    low, high = sorted((start, end))
    starts, ends, prefix = intervals

    def covered_through(position):
        i = bisect_right(starts, position) - 1
        if i < 0:
            return 0
        return prefix[i] + min(position, ends[i]) - starts[i] + 1

    return (covered_through(high) - covered_through(low - 1)) / (high - low + 1)


def _validate_repeatmaskmode(value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('repeatmaskmode must be a finite number between 0 and 1')


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
    _validate_repeatmaskmode(mode)
    if not mode:
        return hsps
    path = paths.get('te')
    subject_te = _load_te_index(path) if path else {}
    if not path:
        print(f"  Subject {paths['fasta']}: no TE annotation; subject-side masking skipped",
              file=sys.stderr)
    return [h for h in hsps
            if not (_te_coverage(query_te, h['qstart'] + query_offset,
                                    h['qend'] + query_offset) >= mode
                    or _te_coverage(subject_te.get(h['chrom']), h['sstart'], h['send']) >= mode)]


def _validate_search(min_hits, window_mult, min_coverage):
    if not isinstance(min_hits, int) or min_hits < 1:
        raise ValueError('min_hits must be a positive integer')
    if not math.isfinite(window_mult) or window_mult <= 0:
        raise ValueError('window_mult must be finite and greater than zero')
    if not 0 <= min_coverage <= 1:
        raise ValueError('min_coverage must be between 0 and 1')


def _terminal_support(hit, query_len, left):
    band = max(1, math.ceil(query_len * 0.10))
    low, high = (1, band) if left else (query_len - band + 1, query_len)
    overlap = max(0, min(high, hit['qend']) - max(low, hit['qstart']) + 1)
    return overlap / band


def _anchor_pairs(region, query_len):
    # Require coverage of at least half of the terminal 10% band, rather
    # than accepting a tiny repeat merely touching the query end.
    left = [h for h in region if _terminal_support(h, query_len, True) >= 0.5]
    right = [h for h in region if _terminal_support(h, query_len, False) >= 0.5]
    for lhit in left:
        for rhit in right:
            if lhit['strand'] != rhit['strand']:
                continue
            if lhit['qstart'] > rhit['qstart'] or lhit['qend'] > rhit['qend']:
                continue
            if lhit['strand'] == '+':
                if lhit['sstart'] > rhit['sstart'] or lhit['send'] > rhit['send']:
                    continue
                low, high = lhit['sstart'], rhit['send']
            else:
                if lhit['sstart'] < rhit['sstart'] or lhit['send'] < rhit['send']:
                    continue
                low, high = rhit['sstart'], lhit['send']
            # Both anchors must themselves be fully inside this interval.
            inside = [h for h in region if low <= h['sstart'] and h['send'] <= high]
            yield inside, lhit, rhit, low, high


def _score_region(region, query_len, status, anchors=None):
    ordered = sorted(region, key=lambda h: h['qstart'])
    coverage = _union(region, 'qstart', 'qend') / query_len
    coll = _collinearity(ordered)
    score = coverage * (0.8 + 0.2 * coll)
    near_start = [h for h in ordered if h['qstart'] <= query_len * 0.10]
    near_end = [h for h in ordered if h['qend'] >= query_len * 0.90]
    if anchors:
        best_start, best_end, low, high = anchors
    else:
        best_start = max(near_start, key=lambda h: h['length']) if near_start else ordered[0]
        best_end = max(near_end, key=lambda h: h['length']) if near_end else ordered[-1]
        low = min(h['sstart'] for h in region)
        high = max(h['send'] for h in region)
    score += max(0, 1.0 - best_start['qstart'] / query_len
                 - (1.0 - best_end['qend'] / query_len)) * 0.15
    score += min(1.0, (best_start['length'] + best_end['length']) / 2000.0) * 0.10
    # Soft penalty, not a span cutoff: retain genuine insertions/expansions.
    span_penalty = 0.05 * math.log2(max(1.0, (high - low + 1) / query_len)) if anchors else 0.0
    score -= span_penalty
    strands = {h['strand'] for h in region}
    weight = sum(h['length'] for h in region)
    return {
        'chrom': region[0]['chrom'], 'start': low, 'end': high,
        'coverage': coverage,
        'identity': sum(h['identity'] * h['length'] for h in region) / weight,
        'score': score,
        'orientation': next(iter(strands)) if len(strands) == 1 else 'mixed',
        'num_hits': len(region), 'subject_span': high - low + 1,
        '_q_low': min(h['qstart'] for h in region),
        '_q_high': max(h['qend'] for h in region),
        'boundary_method': status,
        'anchor_strand': best_start['strand'] if anchors else '.',
        'left_anchor_qstart': best_start['qstart'] if anchors else 0,
        'left_anchor_qend': best_start['qend'] if anchors else 0,
        'left_anchor_sstart': best_start['sstart'] if anchors else 0,
        'left_anchor_send': best_start['send'] if anchors else 0,
        'right_anchor_qstart': best_end['qstart'] if anchors else 0,
        'right_anchor_qend': best_end['qend'] if anchors else 0,
        'right_anchor_sstart': best_end['sstart'] if anchors else 0,
        'right_anchor_send': best_end['send'] if anchors else 0,
        'span_penalty': span_penalty,
    }


def _candidate_rank(candidate):
    return (candidate.get('boundary_method') in ('paired_anchors', 'paired_flank_chains'), candidate['score'])


def _find_candidates(hsps, query_len, min_coverage, window_mult=5, min_hits=2,
                     boundary_mode='chains'):
    _validate_search(min_hits, window_mult, min_coverage)
    if query_len <= 0:
        raise ValueError('Query sequence must not be empty')
    if boundary_mode not in ('chains', 'anchors', 'window'):
        raise ValueError('boundary_mode must be chains, anchors, or window')
    if boundary_mode == 'window':
        return _find_window_candidates(hsps, query_len, min_coverage,
                                       window_mult, min_hits, 'window')
    results = _chain_candidates(hsps, query_len, min_hits, min_coverage)
    if not results:
        # Preserve the previous paired-anchor/window result if no chain pair
        # survives; this is exactly the fallback used in the reviewed trial.
        results = _find_window_candidates(hsps, query_len, min_coverage,
                                         window_mult, min_hits, 'anchors')
        for c in results:
            c['support_status'] = 'legacy_unresolved'
        return results
    for c in results:
        c['support_status'] = ('paired_flank_chains' if c['weak_support'] >= .15
                               else 'weak_flank_support')
    best = results[0]
    alternate = next((c for c in results[1:] if c['chrom'] != best['chrom']
                      or abs(c['start'] - best['start']) > 5000
                      or abs(c['end'] - best['end']) > 5000), None)
    best['score_margin'] = best['score'] - alternate['score'] if alternate else ''
    if alternate and best['score_margin'] < .15:
        best['support_status'] = 'ambiguous_pair'
    return results


def _find_window_candidates(hsps, query_len, min_coverage, window_mult=5, min_hits=2,
                     boundary_mode='chains'):
    _validate_search(min_hits, window_mult, min_coverage)
    if query_len <= 0:
        raise ValueError('Query sequence must not be empty')
    if boundary_mode not in ('chains', 'anchors', 'window'):
        raise ValueError('boundary_mode must be chains, anchors, or window')
    by_chrom = defaultdict(list)
    for h in hsps:
        by_chrom[h['chrom']].append(h)
    candidates = []
    for chrom, group in by_chrom.items():
        if len(group) < min_hits:
            continue
        group.sort(key=lambda h: h['sstart'])
        window = max(1, int(query_len * window_mult))
        step = max(1, window // 2)
        buffer = window // 4
        w_start = group[0]['sstart']
        s_max_all = max(h['send'] for h in group)
        while w_start <= s_max_all:
            seed = [h for h in group if w_start <= h['sstart'] <= w_start + window]
            if len(seed) >= min_hits:
                s_low = min(h['sstart'] for h in seed)
                s_high = max(h['send'] for h in seed)
                region = [h for h in group if s_low - buffer <= h['sstart'] <= s_high + buffer]
                anchored = []
                if boundary_mode == 'anchors':
                    for inside, left, right, low, high in _anchor_pairs(region, query_len):
                        if len(inside) < min_hits:
                            continue
                        c = _score_region(inside, query_len, 'paired_anchors', (left, right, low, high))
                        if c['coverage'] >= min_coverage:
                            c['excluded_hsps'] = len(region) - len(inside)
                            anchored.append(c)
                if anchored:
                    candidates.extend(anchored)
                else:
                    c = _score_region(region, query_len,
                                      'window' if boundary_mode == 'window' else 'window_fallback')
                    c['excluded_hsps'] = 0
                    if c['coverage'] >= min_coverage:
                        candidates.append(c)
            w_start += step
    if boundary_mode == 'anchors':
        # Overlapping windows may rediscover the same anchored interval.
        unique = {}
        for c in candidates:
            key = (c['chrom'], c['start'], c['end'], c['boundary_method'], c['anchor_strand'])
            if key not in unique or _candidate_rank(c) > _candidate_rank(unique[key]):
                unique[key] = c
        candidates = list(unique.values())
    candidates.sort(key=_candidate_rank, reverse=True)
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
    return max(candidates, key=_candidate_rank) if candidates else None


def write_results(all_candidates, best_per_genome, all_path, best_path,
                   query_region=None, genome_order=None):
    diagnostics = ('boundary_method', 'anchor_strand', 'left_anchor_qstart',
                   'left_anchor_qend', 'left_anchor_sstart', 'left_anchor_send',
                   'right_anchor_qstart', 'right_anchor_qend', 'right_anchor_sstart',
                   'right_anchor_send', 'span_penalty', 'excluded_hsps',
                   'support_status', 'left_support', 'right_support',
                   'left_chain_hits', 'right_chain_hits', 'score_margin')
    with open(all_path, 'w') as f:
        f.write("genome\tchrom\tstart\tend\tcoverage\tidentity\tscore\t"
                "orientation\tnum_hits\tsubject_span\tq_low\tq_high\t" + "\t".join(diagnostics) + "\n")
        for c in all_candidates:
            f.write(f"{c['genome']}\t{c['chrom']}\t{c['start']}\t{c['end']}\t"
                    f"{c['coverage']:.4f}\t{c['identity']:.2f}\t{c['score']:.4f}\t"
                    f"{c['orientation']}\t{c['num_hits']}\t{c['subject_span']}\t"
                    f"{c.get('_q_low', 0)}\t{c.get('_q_high', 0)}\t" +
                    "\t".join(str(c.get(key, "")) for key in diagnostics) + "\n")
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
