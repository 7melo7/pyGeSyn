import sys
import subprocess
import tempfile
import os
import platform

import matplotlib
matplotlib.use('Agg')
_os = platform.system()
if _os == 'Windows':
    matplotlib.rcParams['font.family'] = 'Arial'
elif _os == 'Darwin':
    matplotlib.rcParams['font.family'] = 'Helvetica'
else:
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.size'] = 8
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from .genome import extract_sequence
from .gff import parse_gene_structures, parse_te_features
from .config import load_config

COLOR_GENE = '#4878d0'
COLOR_TE = '#ee854a'
COLOR_TRACK = '#999999'
DPI = 200


def draw_dotplot(regions_path, config_path, pair, output_path,
                 min_length=100, min_identity=80.0):
    config = load_config(config_path)
    name_a, name_b = pair.split(',')

    raw = _parse_regions(regions_path)
    region_a = next((r for r in raw if r[0] == name_a), None)
    region_b = next((r for r in raw if r[0] == name_b), None)
    if not region_a or not region_b:
        raise ValueError("Region not found")

    seq_a = extract_sequence(config[name_a]['fasta'],
                             region_a[1], region_a[2], region_a[3])
    seq_b = extract_sequence(config[name_b]['fasta'],
                             region_b[1], region_b[2], region_b[3])

    genes_a = parse_gene_structures(config[name_a]['gff3'],
                                    region_a[1], region_a[2], region_a[3])
    tes_a = parse_te_features(config[name_a]['te'],
                              region_a[1], region_a[2], region_a[3])
    genes_b = parse_gene_structures(config[name_b]['gff3'],
                                    region_b[1], region_b[2], region_b[3])
    tes_b = parse_te_features(config[name_b]['te'],
                              region_b[1], region_b[2], region_b[3])

    hsps = _blast_pair(seq_a, seq_b, min_length, min_identity)

    if len(seq_a) < len(seq_b):
        seq_a, seq_b = seq_b, seq_a
        name_a, name_b = name_b, name_a
        region_a, region_b = region_b, region_a
        genes_a, genes_b = genes_b, genes_a
        tes_a, tes_b = tes_b, tes_a
        for h in hsps:
            h['qstart'], h['sstart'] = h['sstart'], h['qstart']
            h['qend'], h['send'] = h['send'], h['qend']

    len_a = len(seq_a)
    len_b = len(seq_b)

    fig_w = 12
    ratio = len_b / len_a if len_a > 0 else 1
    fig_h = max(5, fig_w * ratio * 0.9 + 2)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)

    gs = fig.add_gridspec(2, 2, width_ratios=[0.06, 0.94],
                          height_ratios=[0.92, 0.08],
                          hspace=0.0, wspace=0.0)

    ax_main = fig.add_subplot(gs[0, 1])
    ax_x = fig.add_subplot(gs[1, 1], sharex=ax_main)
    ax_y = fig.add_subplot(gs[0, 0], sharey=ax_main)

    pos = ax_y.get_position()
    ax_y.set_position([pos.x0, pos.y0, pos.width, pos.height])

    ax_main.set_xlim(0, len_a)
    ax_main.set_ylim(0, len_b)

    feat_range = 1.2
    vscale = feat_range * fig_w / fig_h if fig_h > 0 else 1.2
    ax_x.set_ylim(-feat_range/2, feat_range/2)
    ax_y.set_xlim(-vscale/2, vscale/2)
    for spine in ax_main.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    ax_main.tick_params(labelbottom=False, labelleft=False,
                        bottom=False, left=False)

    for h in hsps:
        ax_main.plot([h['qstart'], h['qend']],
                     [h['sstart'], h['send']],
                     color='black', linewidth=0.6, alpha=0.6,
                     solid_capstyle='butt')
        ax_main.plot([0, h['qstart']], [h['sstart'], h['sstart']],
                     color='red', linewidth=0.2, linestyle=':')
        ax_main.plot([h['qstart'], h['qstart']], [0, h['sstart']],
                     color='red', linewidth=0.2, linestyle=':')
        ax_main.plot([0, h['qend']], [h['send'], h['send']],
                     color='red', linewidth=0.2, linestyle=':')
        ax_main.plot([h['qend'], h['qend']], [0, h['send']],
                     color='red', linewidth=0.2, linestyle=':')

    _draw_track_h(ax_x, genes_a, tes_a, region_a[2])
    _draw_track_v(ax_y, genes_b, tes_b, region_b[2], vscale)

    ax_main.set_xlabel("", fontsize=8, labelpad=4)
    ax_main.text(0.5, 1.004, f"{name_a}  {region_a[1]}:{region_a[2]:,}-{region_a[3]:,}",
                 fontsize=8, ha='center', va='bottom', transform=ax_main.transAxes)
    ax_x.set_xlabel(f"{name_a}  {region_a[1]}:{region_a[2]:,}-{region_a[3]:,}",
                    fontsize=8, labelpad=2)
    ax_main.text(1.002, 0.5, f"{name_b}  {region_b[1]}:{region_b[2]:,}-{region_b[3]:,}",
                 fontsize=8, ha='left', va='center',
                 rotation=-90, transform=ax_main.transAxes)

    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved {output_path}", file=sys.stderr)


def _draw_track_h(ax, genes, tes, region_start):
    ax.axis('off')

    ax.axhline(y=0, color=COLOR_TRACK, linewidth=2.5, zorder=2)

    for t in tes:
        cy = 0.22 if t['strand'] == '+' else -0.22
        _te_h(ax, t, cy, region_start)
    for g in genes:
        cy = 0.22 if g['strand'] == '+' else -0.22
        _gene_h(ax, g, cy, region_start)


def _gene_h(ax, gene, cy, region_start):
    for intron in gene.get('introns', []):
        x1, x2 = intron['start'] - region_start, intron['end'] - region_start + 1
        if x2 > x1:
            ax.plot([x1, x2], [cy, cy], 'k-', linewidth=0.3,
                    solid_capstyle='butt', zorder=4)
    for exon in gene.get('exons', []):
        has_cds = len(exon.get('cds_regions', [])) > 0
        for cs, ce in exon.get('cds_regions', []):
            x1, x2 = cs - region_start, ce - region_start + 1
            if x2 > x1:
                ax.add_patch(patches.Rectangle(
                    (x1, cy - 0.14), x2 - x1, 0.28,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5,
                    clip_on=False))
        for us, ue in exon.get('utr_regions', []):
            x1, x2 = us - region_start, ue - region_start + 1
            if x2 > x1:
                h = 0.28 if not has_cds else 0.14
                ax.add_patch(patches.Rectangle(
                    (x1, cy - h/2), x2 - x1, h,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5,
                    clip_on=False))


def _te_h(ax, te, cy, region_start):
    x1, x2 = te['start'] - region_start, te['end'] - region_start + 1
    if x2 <= x1:
        return
    ax.add_patch(patches.Rectangle(
        (x1, cy - 0.112), x2 - x1, 0.224,
        facecolor=COLOR_TE, linewidth=0, zorder=5, clip_on=False))


def _draw_track_v(ax, genes, tes, region_start, vscale):
    ax.axis('off')
    ax.invert_xaxis()

    ax.axvline(x=0, color=COLOR_TRACK, linewidth=2.5, zorder=2)

    fgap = 0.22 * vscale / 1.2
    fscale = 0.08 / 0.06

    for t in tes:
        cx = fgap if t['strand'] == '+' else -fgap
        _te_v(ax, t, cx, region_start, fscale)
    for g in genes:
        cx = fgap if g['strand'] == '+' else -fgap
        _gene_v(ax, g, cx, region_start, fscale)


def _gene_v(ax, gene, cx, region_start, scale=1.0):
    ch = 0.28 * scale
    uh = 0.14 * scale
    for intron in gene.get('introns', []):
        y1, y2 = intron['start'] - region_start, intron['end'] - region_start + 1
        if y2 > y1:
            ax.plot([cx, cx], [y1, y2], 'k-', linewidth=0.3,
                    solid_capstyle='butt', zorder=4)
    for exon in gene.get('exons', []):
        has_cds = len(exon.get('cds_regions', [])) > 0
        for cs, ce in exon.get('cds_regions', []):
            y1, y2 = cs - region_start, ce - region_start + 1
            if y2 > y1:
                ax.add_patch(patches.Rectangle(
                    (cx - ch/2, y1), ch, y2 - y1,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5,
                    clip_on=False))
        for us, ue in exon.get('utr_regions', []):
            y1, y2 = us - region_start, ue - region_start + 1
            if y2 > y1:
                h = ch if not has_cds else uh
                ax.add_patch(patches.Rectangle(
                    (cx - h/2, y1), h, y2 - y1,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5,
                    clip_on=False))


def _te_v(ax, te, cx, region_start, scale=1.0):
    th = 0.224 * scale
    y1, y2 = te['start'] - region_start, te['end'] - region_start + 1
    if y2 <= y1:
        return
    ax.add_patch(patches.Rectangle(
        (cx - th/2, y1), th, y2 - y1,
        facecolor=COLOR_TE, linewidth=0, zorder=5, clip_on=False))


def _blast_pair(seq_a, seq_b, min_length, min_identity):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fa:
        fa.write(f'>a\n{seq_a}\n')
        fa_path = fa.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fb:
        fb.write(f'>b\n{seq_b}\n')
        fb_path = fb.name
    try:
        r = subprocess.run(
            ['blastn', '-query', fa_path, '-subject', fb_path,
             '-outfmt', '6 qstart qend sstart send length pident',
             '-evalue', '1e-5'],
            capture_output=True, text=True)
        hsps = []
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            p = line.split('\t')
            if len(p) < 6:
                continue
            ln, pid = int(p[4]), float(p[5])
            if ln < min_length or pid < min_identity:
                continue
            qs, qe = int(p[0]), int(p[1])
            ss, se = int(p[2]), int(p[3])
            if qs > qe: qs, qe = qe, qs
            if ss > se: ss, se = se, ss
            hsps.append({'qstart': qs, 'qend': qe, 'sstart': ss, 'send': se})
        return hsps
    finally:
        os.unlink(fa_path)
        os.unlink(fb_path)


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
