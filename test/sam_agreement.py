#!/usr/bin/env python3
"""Per-read agreement between two SAM files, primary records only.

    python3 test/sam_agreement.py star_Aligned.out.sam rustar_Aligned.out.sam

Reports how many mates land at the same chromosome, position and CIGAR, and
how many carry the same NH. Ties broken differently by the two aligners show
up as position differences with equal NH, which is why the two figures are
printed apart: a drop in NH agreement is a different problem from a drop in
position agreement.
"""
import sys
def load(path):
    out={}
    for line in open(path):
        if line.startswith("@"): continue
        f=line.rstrip("\n").split("\t")
        flag=int(f[1])
        if flag & 0x100 or flag & 0x800: continue
        mate = 2 if flag & 0x80 else 1
        key=(f[0], mate)
        nh=0
        for t in f[11:]:
            if t.startswith("NH:i:"): nh=int(t[5:])
        out[key]=(f[2], f[3], f[5], nh, flag & 0x4)
    return out
a,b=load(sys.argv[1]),load(sys.argv[2])
keys=set(a)|set(b)
same=pos_same=nh_same=0
only_a=only_b=0
nh_hist_a={}; nh_hist_b={}
for k in keys:
    x,y=a.get(k),b.get(k)
    if x is None: only_b+=1; continue
    if y is None: only_a+=1; continue
    if x[4]==0: nh_hist_a[x[3]]=nh_hist_a.get(x[3],0)+1
    if y[4]==0: nh_hist_b[y[3]]=nh_hist_b.get(y[3],0)+1
    if x[:3]==y[:3]: pos_same+=1
    if x[:4]==y[:4]: same+=1
    if x[3]==y[3]: nh_same+=1
n=len(keys)
print(f"records compared: {n}   only in A: {only_a}   only in B: {only_b}")
print(f"same chr/pos/CIGAR : {pos_same} ({pos_same/n:.4%})")
print(f"same incl. NH      : {same} ({same/n:.4%})")
print(f"same NH            : {nh_same} ({nh_same/n:.4%})")
print("max NH  A:", max(nh_hist_a or {0:0}), " B:", max(nh_hist_b or {0:0}))
