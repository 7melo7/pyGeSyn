"""Terminal-seeded flank chaining validated on the GM4 310-genome review.

Scores and support flags are heuristics, not probabilities of orthology.
Multiplicity affects ranking only; it never removes an HSP.
"""
import collections
import math


def _union(items, sk, ek):
    intervals = sorted((h[sk], h[ek]) for h in items)
    total, end = 0, -1
    for start, stop in intervals:
        total += max(0, stop - max(start, end + 1) + 1)
        end = max(end, stop)
    return total


def flank_chains(hits, left, query_len):
    L = query_len
    # Work outwards from a terminal HSP into its flank, independently per strand.
    band=max(1,int(.25*L))
    eligible=[h for h in hits if (h['qstart']<=band if left else h['qend']>=L-band+1)]
    def coords(h):
        qs,qe=h['qstart'],h['qend']
        ss,se=(h['sstart'],h['send']) if h['strand']=='+' else (-h['send'],-h['sstart'])
        if not left: qs,qe,ss,se=L-qe+1,L-qs+1,-se,-ss
        return qs,qe,ss,se
    eligible.sort(key=lambda h:coords(h)[0])
    coords_list=[coords(h) for h in eligible]
    chains=[]
    for a,anchor in enumerate(eligible):
        aq,ae,_,_=coords_list[a]
        if aq>max(1000,int(.02*L)): continue
        # Query-distance bins approximate multi-locus support, without counting
        # overlapping HSPs at the same subject locus as independent copies.
        def contribution(h,start,end):
            return max(0,min(end,band)-max(start,1)+1)*h['identity']/100/h['_copies']**.5
        scores={a:contribution(anchor,aq,ae)}; paths={a:[a]}
        for j in range(a+1,len(eligible)):
            h=eligible[j]; q0,q1,s0,s1=coords_list[j]
            if h['strand']!=anchor['strand']: continue
            best=None
            for i,value in list(scores.items()):
                iq,ie,iss,ise=coords_list[i]
                if q0<=iq or q1<=ie or s0<=iss or s1<=ise: continue
                qgap=q0-ie-1; sgap=s0-ise-1
                if qgap>band*.6 or sgap>band: continue
                if min(qgap,sgap)<-min(ie-iq+1,q1-q0+1)*.5: continue
                gain=contribution(h,max(q0,ie+1),q1)
                penalty=.15*abs(sgap-qgap)+.02*max(0,qgap,sgap)
                value2=value+gain-penalty
                if best is None or value2>best[0]:best=(value2,i)
            if best is not None and best[0]>0:
                scores[j]=best[0];paths[j]=paths[best[1]]+[j]
        last=max(scores,key=scores.get); chain=[eligible[i] for i in paths[last]]
        # Terminal proximity matters independently of flank coverage.
        terminal=max(0,1-(aq-1)/max(1000,.02*L))
        support=min(1,scores[last]/(L*.08))*terminal
        chains.append(dict(anchor=anchor,support=support,chain=chain,
                           bases=_union(chain,'qstart','qend'),terminal_distance=aq-1))
    return chains

def candidates(hits, query_len, min_hits=2, min_coverage=0):
    L = query_len
    hits = [dict(h) for h in hits]
    # Estimate repeated loci for each query segment. Used for scoring only.
    for h in hits:
        mid=(h['qstart']+h['qend'])//2
        loci={(x['chrom'],x['sstart']//10000) for x in hits
              if x['qstart']<=mid<=x['qend'] and x['identity']>=h['identity']-3}
        h['_copies']=max(1,len(loci))
    groups=collections.defaultdict(list)
    for h in hits:groups[h['chrom']].append(h)
    result={}
    for chrom,group in groups.items():
        for left in flank_chains(group,True,L):
            for right in flank_chains(group,False,L):
                a,b=left['anchor'],right['anchor']
                if a['strand']!=b['strand'] or a['qstart']>b['qstart'] or a['qend']>b['qend']:continue
                if a['strand']=='+':
                    if a['sstart']>b['sstart'] or a['send']>b['send']:continue
                    low,high=a['sstart'],b['send']
                else:
                    if a['sstart']<b['sstart'] or a['send']<b['send']:continue
                    low,high=b['sstart'],a['send']
                if any(h['sstart']<low or h['send']>high for h in left['chain']+right['chain']):continue
                inside=[h for h in group if low<=h['sstart'] and h['send']<=high]
                if len(inside)<min_hits:continue
                cov=_union(inside,'qstart','qend')/L
                if cov < min_coverage: continue
                weak=min(left['support'],right['support'])
                score=2*weak+.5*(left['support']+right['support'])+.15*cov
                score-=.1*math.log2(max(1,(high-low+1)/L))
                c=dict(chrom=chrom,start=low,end=high,strand=a['strand'],score=score,
                       coverage=cov,left_support=left['support'],right_support=right['support'],
                       left_chain_hits=len(left['chain']),right_chain_hits=len(right['chain']),
                       left_anchor_qstart=a['qstart'],left_anchor_qend=a['qend'],
                       right_anchor_qstart=b['qstart'],right_anchor_qend=b['qend'],
                       num_hits=len(inside),weak_support=weak,
                       left_chain=left['chain'],right_chain=right['chain'])
                weight = sum(h['length'] for h in inside)
                strands = {h['strand'] for h in inside}
                c.update(boundary_method='paired_flank_chains', anchor_strand=a['strand'],
                         identity=sum(h['identity']*h['length'] for h in inside)/weight,
                         orientation=next(iter(strands)) if len(strands)==1 else 'mixed',
                         subject_span=high-low+1, _q_low=min(h['qstart'] for h in inside),
                         _q_high=max(h['qend'] for h in inside),
                         left_anchor_sstart=a['sstart'], left_anchor_send=a['send'],
                         right_anchor_sstart=b['sstart'], right_anchor_send=b['send'],
                         span_penalty=.1*math.log2(max(1,(high-low+1)/L)),
                         excluded_hsps=len(group)-len(inside))
                key=(chrom,low,high,a['strand'])
                if key not in result or score>result[key]['score']:result[key]=c
    return sorted(result.values(),key=lambda c:c['score'],reverse=True)

