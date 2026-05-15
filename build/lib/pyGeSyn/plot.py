import sys
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
matplotlib.rcParams['font.size'] = 6
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath


FEATURE_HEIGHT = 20
UTR_HEIGHT = 10
TE_HEIGHT = 16
FEATURE_GAP = 15
CONN_PAD = FEATURE_GAP - FEATURE_HEIGHT / 2

COLOR_GENE = '#4878d0'
COLOR_TE = '#ee854a'
COLOR_CONN_FWD = '#D9D9D9'
COLOR_CONN_INV = '#FFB3B3'
ALPHA_CONN = 0.65

A4_W_INCH = 8.27
A4_H_INCH = 11.69
DPI = 200


def draw_synteny(regions, all_features, blast_results, output_path,
                 conn_ratio=1.5):
    n = len(regions)
    region_lengths = [r[3] - r[2] + 1 for r in regions]
    max_len = max(region_lengths)

    label_frac = 0.10
    label_pad_px = max_len * 0.015

    xlim_left = -max_len * label_frac - label_pad_px
    xlim_right = max_len * 1.02
    x_range = xlim_right - xlim_left

    track_x_range = max_len
    track_frac_of_axes = track_x_range / x_range

    target_track_inch = A4_W_INCH * 0.78
    axes_w_inch = target_track_inch / track_frac_of_axes
    left_inch = A4_W_INCH * 0.10
    right_inch = A4_W_INCH - left_inch - axes_w_inch

    feature_extent = FEATURE_GAP + FEATURE_HEIGHT / 2
    conn_zone = conn_ratio * feature_extent
    track_spacing = 2 * feature_extent + conn_zone
    axes_h_px = track_spacing * (n - 1) + 2 * feature_extent + 60
    axes_h_inch = axes_h_px / DPI

    top_inch = 0.15
    bottom_margin = 0.60
    axes_bottom_inch = max(bottom_margin, A4_H_INCH - top_inch - axes_h_inch)
    usable_h = A4_H_INCH - top_inch - bottom_margin

    if axes_h_inch > usable_h:
        scale = usable_h / axes_h_inch
        FEATURE_HEIGHT_s = FEATURE_HEIGHT * scale
        UTR_HEIGHT_s = UTR_HEIGHT * scale
        TE_HEIGHT_s = TE_HEIGHT * scale
        FEATURE_GAP_s = FEATURE_GAP * scale
        CONN_PAD_s = CONN_PAD * scale
        feature_extent_s = FEATURE_GAP_s + FEATURE_HEIGHT_s / 2
        conn_zone_s = conn_ratio * feature_extent_s
        track_spacing_s = 2 * feature_extent_s + conn_zone_s
        axes_h_px_s = track_spacing_s * (n - 1) + 2 * feature_extent_s + 60
        axes_h_inch = axes_h_px_s / DPI
    else:
        FEATURE_HEIGHT_s = FEATURE_HEIGHT
        UTR_HEIGHT_s = UTR_HEIGHT
        TE_HEIGHT_s = TE_HEIGHT
        FEATURE_GAP_s = FEATURE_GAP
        CONN_PAD_s = CONN_PAD
        feature_extent_s = feature_extent
        conn_zone_s = conn_zone
        track_spacing_s = track_spacing
        axes_h_px_s = axes_h_px

    fig = plt.figure(figsize=(A4_W_INCH, A4_H_INCH), dpi=DPI)
    ax = fig.add_axes([
        left_inch / A4_W_INCH,
        axes_bottom_inch / A4_H_INCH,
        axes_w_inch / A4_W_INCH,
        axes_h_inch / A4_H_INCH,
    ])

    y_centers = [
        axes_h_px_s - feature_extent_s - i * track_spacing_s
        for i in range(n)
    ]

    for i, region in enumerate(regions):
        name, chrom, start, end, seq = region
        rlen = region_lengths[i]
        y = y_centers[i]
        ax.add_patch(patches.Rectangle(
            (0, y - 1.5), rlen - 1, 3,
            facecolor='#8b8b8b', edgecolor='#8b8b8b',
            linewidth=0.6, zorder=3))
        _draw_features(ax, all_features[i], y, start,
                       FEATURE_HEIGHT_s, UTR_HEIGHT_s, TE_HEIGHT_s,
                       FEATURE_GAP_s)
        _draw_track_label(ax, name, chrom, start, end, rlen, y)

    for pair_idx, hits in enumerate(blast_results):
        if not hits:
            continue
        y_top = y_centers[pair_idx]
        y_bottom = y_centers[pair_idx + 1]
        _draw_blast_connections(ax, hits, y_top, y_bottom,
                                feature_extent_s, CONN_PAD_s)

    ax.set_xlim(xlim_left, xlim_right)
    ax.set_ylim(-20, axes_h_px_s + 20)
    ax.set_yticks([])
    for spine in ['top', 'right', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)

    _add_legend(ax)

    plt.savefig(output_path, dpi=DPI)
    plt.close()
    print(f"  Plot saved to {output_path}", file=sys.stderr)


def _draw_features(ax, features, track_y, region_start,
                   fh, uh, th, fg):
    for feat in features:
        category = feat.get('category', 'gene')
        strand = feat['strand']
        if strand == '+':
            cy = track_y + fg
        else:
            cy = track_y - fg
        if category == 'gene':
            _draw_gene(ax, feat, cy, region_start, fh, uh)
        elif category == 'te':
            _draw_te(ax, feat, cy, region_start, th)


def _draw_gene(ax, gene, center_y, region_start, fh, uh):
    for intron in gene.get('introns', []):
        x1 = intron['start'] - region_start
        x2 = intron['end'] - region_start + 1
        if x2 > x1:
            ax.plot([x1, x2], [center_y, center_y],
                    'k-', linewidth=0.5, solid_capstyle='butt', zorder=4)

    for exon in gene.get('exons', []):
        has_cds = len(exon.get('cds_regions', [])) > 0

        for cs, ce in exon.get('cds_regions', []):
            x1 = cs - region_start
            x2 = ce - region_start + 1
            w = x2 - x1
            if w > 0:
                ax.add_patch(patches.Rectangle(
                    (x1, center_y - fh / 2), w, fh,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5))

        for us, ue in exon.get('utr_regions', []):
            x1 = us - region_start
            x2 = ue - region_start + 1
            w = x2 - x1
            if w > 0:
                h = fh if not has_cds else uh
                ax.add_patch(patches.Rectangle(
                    (x1, center_y - h / 2), w, h,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5))

        if not has_cds and not exon.get('utr_regions') and not exon.get('cds_regions'):
            x1 = exon['exon_start'] - region_start
            x2 = exon['exon_end'] - region_start + 1
            w = x2 - x1
            if w > 0:
                ax.add_patch(patches.Rectangle(
                    (x1, center_y - fh / 2), w, fh,
                    facecolor=COLOR_GENE, linewidth=0, zorder=5))


def _draw_te(ax, te, center_y, region_start, th):
    x1 = te['start'] - region_start
    x2 = te['end'] - region_start + 1
    w = x2 - x1
    if w <= 0:
        return
    ax.add_patch(patches.Rectangle(
        (x1, center_y - th / 2), w, th,
        facecolor=COLOR_TE, linewidth=0, zorder=5))


def _draw_track_label(ax, name, chrom, start, end, rlen, y):
    label = f"{name}  {chrom}:{start:,}-{end:,}"
    ax.text(-rlen * 0.015, y, label, ha='right', va='center',
            fontweight='bold')


def _draw_blast_connections(ax, hits, y_top, y_bottom,
                             feature_extent, conn_pad):
    y_conn_upper = y_top - feature_extent - conn_pad
    y_conn_lower = y_bottom + feature_extent + conn_pad
    mid_y = (y_conn_upper + y_conn_lower) / 2

    for hit in hits:
        qstart = hit['qstart']
        qend = hit['qend']
        sstart = hit['sstart']
        send = hit['send']

        q_left = qstart - 1
        q_right = qend

        if sstart <= send:
            s_left = sstart - 1
            s_right = send
            color = COLOR_CONN_FWD
        else:
            s_left = send
            s_right = sstart - 1
            color = COLOR_CONN_INV

        verts = [
            (q_left,  y_conn_upper),
            (q_left,  mid_y),
            (s_left,  mid_y),
            (s_left,  y_conn_lower),
            (s_right, y_conn_lower),
            (s_right, mid_y),
            (q_right, mid_y),
            (q_right, y_conn_upper),
            (q_left,  y_conn_upper),
        ]
        codes = [
            mpath.Path.MOVETO,
            mpath.Path.CURVE4, mpath.Path.CURVE4, mpath.Path.CURVE4,
            mpath.Path.LINETO,
            mpath.Path.CURVE4, mpath.Path.CURVE4, mpath.Path.CURVE4,
            mpath.Path.CLOSEPOLY,
        ]
        ax.add_patch(patches.PathPatch(
            mpath.Path(verts, codes),
            facecolor=color, edgecolor='black',
            linewidth=0.3, alpha=ALPHA_CONN, zorder=1))


def _add_legend(ax):
    from matplotlib.lines import Line2D
    handles = [
        patches.Patch(facecolor=COLOR_GENE, label='Gene'),
        patches.Patch(facecolor=COLOR_TE, label='TE'),
        Line2D([0], [0], color='black', linewidth=0.5, label='Intron'),
        patches.Patch(facecolor=COLOR_CONN_FWD, label='Collinear'),
        patches.Patch(facecolor=COLOR_CONN_INV, label='Inverted'),
    ]
    ax.legend(
        handles=handles, loc='upper center',
        bbox_to_anchor=(0.5, -0.06), ncol=5,
        frameon=False, fontsize=6)
