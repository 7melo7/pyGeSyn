def _parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _clip(val, lo, hi):
    return max(lo, min(val, hi))


def parse_gene_structures(gff_path, chrom, start, end):
    features = []
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            seqid = parts[0]
            ftype = parts[2]
            try:
                fstart = int(parts[3])
                fend = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            attrs = _parse_attributes(parts[8])

            if seqid != chrom:
                continue
            if fend < start or fstart > end:
                continue

            features.append({
                'type': ftype,
                'start': fstart,
                'end': fend,
                'strand': strand,
                'id': attrs.get('ID', ''),
                'parent': attrs.get('Parent', ''),
            })

    gene_dict = {}
    mrna_list = []
    exon_by_parent = {}
    cds_by_parent = {}

    for feat in features:
        ftype = feat['type']
        if ftype == 'gene':
            gene_dict[feat['id']] = feat
        elif ftype in ('mRNA', 'transcript'):
            mrna_list.append(feat)
        elif ftype == 'exon':
            p = feat['parent']
            exon_by_parent.setdefault(p, []).append(feat)
        elif ftype == 'CDS':
            p = feat['parent']
            cds_by_parent.setdefault(p, []).append(feat)

    gene_structures = []

    for gene_id, gene in gene_dict.items():
        children = [m for m in mrna_list if m['parent'] == gene_id]
        children.sort(key=lambda m: m['start'])
        if not children:
            continue

        mrna = children[0]
        mid = mrna['id']
        raw_exons = sorted(exon_by_parent.get(mid, []), key=lambda e: e['start'])
        raw_cds_list = sorted(cds_by_parent.get(mid, []), key=lambda c: c['start'])

        exons = []
        for ex in raw_exons:
            e_start = _clip(ex['start'], start, end)
            e_end = _clip(ex['end'], start, end)
            if e_start > e_end:
                continue

            cds_regions = []
            for cd in raw_cds_list:
                if cd['end'] < e_start or cd['start'] > e_end:
                    continue
                cs = _clip(cd['start'], e_start, e_end)
                ce = _clip(cd['end'], e_start, e_end)
                if cs <= ce:
                    cds_regions.append((cs, ce))

            utr_regions = []
            pos = e_start
            for cs, ce in cds_regions:
                if pos < cs:
                    utr_regions.append((pos, cs - 1))
                pos = ce + 1
            if pos <= e_end:
                utr_regions.append((pos, e_end))

            exons.append({
                'exon_start': e_start,
                'exon_end': e_end,
                'cds_regions': cds_regions,
                'utr_regions': utr_regions,
            })

        introns = []
        for i in range(len(exons) - 1):
            intron_start = exons[i]['exon_end'] + 1
            intron_end = exons[i + 1]['exon_start'] - 1
            if intron_start <= intron_end:
                introns.append({
                    'start': _clip(intron_start, start, end),
                    'end': _clip(intron_end, start, end),
                })

        gene_structures.append({
            'category': 'gene',
            'strand': gene['strand'],
            'gene_start': _clip(gene['start'], start, end),
            'gene_end': _clip(gene['end'], start, end),
            'gene_id': gene_id,
            'exons': exons,
            'introns': introns,
        })

    gene_structures.sort(key=lambda g: g['gene_start'])
    return gene_structures


def parse_te_features(gff_path, chrom, start, end):
    features = []
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            seqid = parts[0]
            ftype = parts[2]
            try:
                fstart = int(parts[3])
                fend = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            attrs = _parse_attributes(parts[8])

            if seqid != chrom:
                continue
            if fend < start or fstart > end:
                continue

            clip_start = max(fstart, start)
            clip_end = min(fend, end)

            features.append({
                'category': 'te',
                'strand': strand,
                'start': clip_start,
                'end': clip_end,
                'id': attrs.get('ID', attrs.get('Name', '')),
                'type': ftype,
            })

    features.sort(key=lambda x: x['start'])
    return features
