#!/usr/bin/env python3
import re, sys
import pandas as pd
sam_path, reference_label, out_path = sys.argv[1:4]
rows=[]
def ref_consumed(cigar):
    return sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MDN=X')
with open(sam_path) as h:
    for line in h:
        if line.startswith('@'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<11: continue
        qname, flag_s, rname, pos_s, mapq_s, cigar = f[:6]
        flag=int(flag_s)
        if flag & 0x4: continue
        tags={}
        for tag in f[11:]:
            b=tag.split(':',2)
            if len(b)==3: tags[b[0]]=b[2]
        nm=int(tags.get('NM',-1))
        if nm!=0: continue
        start0=int(pos_s)-1
        end0=start0+ref_consumed(cigar)
        rows.append({'spacer_id':qname,'reference':reference_label,'contig':rname,
                     'start0':start0,'end0':end0,'strand':'-' if flag & 0x10 else '+',
                     'CIGAR':cigar,'MAPQ':int(mapq_s),'NM':nm,
                     'secondary':bool(flag & 0x100),'supplementary':bool(flag & 0x800)})
cols=['spacer_id','reference','contig','start0','end0','strand','CIGAR','MAPQ','NM','secondary','supplementary']
pd.DataFrame(rows,columns=cols).to_csv(out_path,sep='\t',index=False,compression='gzip')
print(f'{reference_label}: {len(rows):,} exact 19-nt alignment records')
