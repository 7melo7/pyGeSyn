"""Gene-level links from longest-CDS transcripts, independent of genomic BLAST."""
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

from .config import load_config, load_gene_colors
from .genome import extract_sequence
from .gff import _parse_attributes, parse_gene_structures, parse_te_features
from .plot import draw_synteny
from .workdir import create_workdir


def union(intervals):
    merged = []
    for start, end in sorted((a, b) for a, b in intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    return merged


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn',
                                      'TGCAYRKMVHDBNtgcayrkmvhdbn'))[::-1]


def representatives(paths, region, workdir, track):
    name, chrom, start, end = region
    records = []
    with open(paths['gff3']) as stream:
        for line in stream:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[0] != chrom:
                continue
            a, b = int(fields[3]), int(fields[4])
            if a < 1 or b < a:
                raise ValueError('Invalid GFF coordinates')
            records.append((fields, _parse_attributes(fields[8])))
    genes = {attr['ID']: fields for fields, attr in records
             if fields[2] == 'gene' and 'ID' in attr
             and int(fields[3]) <= end and int(fields[4]) >= start}
    transcripts, cds, exons = {}, defaultdict(list), defaultdict(list)
    for fields, attr in records:
        kind = fields[2]
        if kind in ('mRNA', 'transcript') and 'ID' in attr:
            for parent in attr.get('Parent', '').split(','):
                if parent in genes:
                    transcripts.setdefault(parent, []).append((attr['ID'], fields))
        elif kind in ('CDS', 'exon'):
            for parent in attr.get('Parent', '').split(','):
                (cds if kind == 'CDS' else exons)[parent].append((int(fields[3]), int(fields[4])))
    # Read once per region, including full selected gene spans beyond the display.
    lo = min([start] + [int(g[3]) for g in genes.values()])
    hi = max([end] + [int(g[4]) for g in genes.values()])
    sequence = extract_sequence(paths['fasta'], chrom, lo, hi)
    if len(sequence) != hi-lo+1:
        raise ValueError(f'{name}/{chrom}: annotation or region extends beyond FASTA')
    colors = load_gene_colors(paths['gene_colors']) if 'gene_colors' in paths else {}
    features = parse_gene_structures(paths['gff3'], chrom, start, end, colors)
    structures = {g['gene_id']: g for g in features}
    result, manifest = [], []
    for gid, gene in sorted(genes.items(), key=lambda x: (int(x[1][3]), x[0])):
        choices = [(mid, fields, union(cds[mid])) for mid, fields in transcripts.get(gid, []) if cds[mid]]
        if not choices:
            manifest.append(dict(track=track, genome=name, gene_id=gid, transcript_id='',
                                 cds_length=0, status='no_cds', boundary_clipped=False))
            continue
        mid, tx, intervals = min(choices, key=lambda x: (-sum(b-a+1 for a,b in x[2]), x[0]))
        if tx[6] not in ('+', '-'):
            raise ValueError(f'{gid}/{mid}: CDS strand must be + or -')
        if any(a < lo or b > hi for a,b in intervals):
            raise ValueError(f'{gid}/{mid}: CDS outside annotated gene bounds')
        seq = ''.join(sequence[a-lo:b-lo+1] for a,b in intervals)
        if tx[6] == '-':
            seq = reverse_complement(seq)
        # Build the displayed structure for the same representative transcript.
        raw_exons = union(exons[mid] + intervals)
        displayed = []
        for a,b in raw_exons:
            a,b = max(a,start),min(b,end)
            if a>b: continue
            coding = [(max(a,c),min(b,d)) for c,d in intervals if c<=b and d>=a]
            utr, pos = [], a
            for c,d in coding:
                if pos<c: utr.append((pos,c-1))
                pos=d+1
            if pos<=b: utr.append((pos,b))
            displayed.append(dict(exon_start=a,exon_end=b,cds_regions=coding,utr_regions=utr))
        structure = dict(category='gene',gene_id=gid,strand=tx[6],
                         gene_start=max(start,int(gene[3])),gene_end=min(end,int(gene[4])),
                         exons=displayed,introns=[dict(start=x['exon_end']+1,end=y['exon_start']-1)
                          for x,y in zip(displayed,displayed[1:]) if y['exon_start']>x['exon_end']+1])
        # Existing parser resolves GeneID/any-isoform colors, including conflicts.
        if gid in structures and 'color' in structures[gid]:
            structure['color']=structures[gid]['color']
        structures[gid]=structure
        key=f't{track}_g{len(result)}'
        clipped=int(gene[3])<start or int(gene[4])>end
        result.append(dict(id=key,gene_id=gid,transcript_id=mid,seq=seq,strand=tx[6],
                           low=intervals[0][0],high=intervals[-1][1],clipped=clipped))
        manifest.append(dict(track=track,genome=name,gene_id=gid,transcript_id=mid,
                             cds_length=len(seq),status='selected',boundary_clipped=clipped))
    fasta=workdir/f'track_{track}.cds.fa'
    fasta.write_text(''.join(f">{r['id']}\n{r['seq']}\n" for r in result))
    return result, list(structures.values()), manifest, fasta


def best_chain(hsps):
    """Max-bit-score compatible chain; overlapping HSPs are alternatives."""
    best = []
    best_score = -1
    for strand in ('+', '-'):
        hits = sorted((h for h in hsps if h['strand']==strand), key=lambda h:(h['qs'],h['qe']))
        scores, paths = [], []
        for j,h in enumerate(hits):
            score, path = h['bits'], [h]
            for i,p in enumerate(hits[:j]):
                ordered = p['se'] < h['ss'] if strand=='+' else p['ss'] > h['se']
                if p['qe'] < h['qs'] and ordered and scores[i]+h['bits']>score:
                    score,path=scores[i]+h['bits'],paths[i]+[h]
            scores.append(score);paths.append(path)
            if score>best_score: best_score,best=score,path
    return best


def compare(a, b, fa, fb, workdir, pair_index, min_identity, min_coverage):
    raw=workdir/f'pair_{pair_index}.blast.tsv'
    if not a or not b:
        raw.write_text('');return []
    command=['blastn','-query',str(fa),'-subject',str(fb),'-evalue','1e-5',
             '-max_target_seqs',str(max(500,len(b))),'-outfmt',
             '6 qseqid sseqid qstart qend sstart send length nident bitscore']
    completed=subprocess.run(command,capture_output=True,text=True)
    if completed.returncode:
        raise RuntimeError(f'blastn failed: {completed.stderr}')
    raw.write_text(completed.stdout)
    grouped=defaultdict(list)
    for line in completed.stdout.splitlines():
        q,s,qs,qe,ss,se,length,nident,bits=line.split()
        qs,qe,ss,se=map(int,(qs,qe,ss,se))
        if qs>qe: qs,qe,ss,se=qe,qs,se,ss
        grouped[q,s].append(dict(qs=qs,qe=qe,ss=min(ss,se),se=max(ss,se),
                                strand='+' if ss<=se else '-',length=int(length),
                                nident=int(nident),bits=float(bits)))
    qa,sb={x['id']:x for x in a},{x['id']:x for x in b}
    results=[]
    for (q,s),hsps in sorted(grouped.items()):
        chain=best_chain(hsps)
        identity=100*sum(h['nident'] for h in chain)/sum(h['length'] for h in chain)
        qc=sum(h['qe']-h['qs']+1 for h in chain)/len(qa[q]['seq'])
        sc=sum(h['se']-h['ss']+1 for h in chain)/len(sb[s]['seq'])
        results.append(dict(query=q,subject=s,identity=identity,query_coverage=qc,
                            subject_coverage=sc,hsp_count=len(chain),cds_strand=chain[0]['strand'],
                            passed=identity>=min_identity and max(qc,sc)>=min_coverage))
    return results


def run(config_path, raw_regions, output, min_identity=90, min_coverage=.8,
        conn_ratio=1.5, pairs_output=None):
    if not math.isfinite(min_identity) or not 0<=min_identity<=100:
        raise ValueError('min-identity must be between 0 and 100')
    if not math.isfinite(min_coverage) or not 0<=min_coverage<=1:
        raise ValueError('min-cds-coverage must be between 0 and 1')
    if not raw_regions: raise ValueError('At least one region is required')
    config=load_config(config_path)
    workdir=create_workdir()
    regions,features,reps,fastas,manifest=[],[],[],[],[]
    for i,r in enumerate(raw_regions):
        name,chrom,start,end=r
        if start<1 or end<start: raise ValueError('Invalid region coordinates')
        paths=config[name]
        if 'gff3' not in paths: raise ValueError(f'{name}: plot_cds requires gff3')
        selected,genes,rows,fa=representatives(paths,r,workdir,i)
        tes=parse_te_features(paths['te'],chrom,start,end) if 'te' in paths else []
        regions.append((*r,''));features.append(genes+tes)
        reps.append(selected);fastas.append(fa);manifest.extend(rows)
    all_pairs,connections=[],[]
    for i in range(len(regions)-1):
        matches=compare(reps[i],reps[i+1],fastas[i],fastas[i+1],workdir,i,min_identity,min_coverage)
        qa,sb={x['id']:x for x in reps[i]},{x['id']:x for x in reps[i+1]}
        links=[]
        for m in matches:
            a,b=qa[m['query']],sb[m['subject']]
            orientation='+' if (a['strand']==b['strand'])==(m['cds_strand']=='+') else '-'
            all_pairs.append(dict(pair_index=i,query_genome=regions[i][0],subject_genome=regions[i+1][0],
                query_gene=a['gene_id'],subject_gene=b['gene_id'],query_transcript=a['transcript_id'],
                subject_transcript=b['transcript_id'],query_cds_length=len(a['seq']),subject_cds_length=len(b['seq']),
                query_boundary_clipped=a['clipped'],subject_boundary_clipped=b['clipped'],orientation=orientation,**m))
            if not m['passed']: continue
            ra,rb=regions[i],regions[i+1]
            ql,qr=max(a['low'],ra[2])-ra[2]+1,min(a['high'],ra[3])-ra[2]+1
            sl,sr=max(b['low'],rb[2])-rb[2]+1,min(b['high'],rb[3])-rb[2]+1
            if ql>qr or sl>sr: continue  # CDS entirely outside visible region.
            links.append(dict(qstart=ql,qend=qr,sstart=sl if orientation=='+' else sr,
                              send=sr if orientation=='+' else sl))
        connections.append(links)
    fields=['pair_index','query_genome','subject_genome','query_gene','subject_gene','query_transcript',
            'subject_transcript','query_cds_length','subject_cds_length','identity','query_coverage',
            'subject_coverage','hsp_count','orientation','passed','query_boundary_clipped','subject_boundary_clipped',
            'query','subject','cds_strand']
    destination=Path(pairs_output) if pairs_output else Path(output).with_suffix('.gene_pairs.tsv')
    for path,rows,columns in [(destination,all_pairs,fields),
                              (workdir/'representatives.tsv',manifest,
                               ['track','genome','gene_id','transcript_id','cds_length','status','boundary_clipped'])]:
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,columns,delimiter='\t');writer.writeheader();writer.writerows(rows)
    draw_synteny(regions,features,connections,output,conn_ratio=conn_ratio)
    print(f'Gene pairs: {destination}\nIntermediate files: {workdir.resolve()}')
    return destination,workdir
