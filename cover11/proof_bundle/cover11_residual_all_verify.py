#!/usr/bin/env python3
"""Independent replay verifier for all 53 noncandidate residual topology certificates."""
from __future__ import annotations
from fractions import Fraction as F
import argparse,glob,json,os,re,sys,time
sys.path.insert(0,'/mnt/data')
import cover11_residual_exact_worker as W
import cover11_candidate_exact_core as core
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int,default=53);ap.add_argument('--out');a=ap.parse_args()
# Corrected T+ metric residual set.
M=json.load(open('/mnt/data/cover11_metric_Tplus.json'))
metric={}
for key,b in M['cases'].items():
 B,I=map(int,key.split(','))
 for z in b['residual']:metric[(B,I,int(z['idx']))]=tuple(tuple(map(int,f)) for f in z['faces'])
candidate=(9,2,2547);wheel=(10,1,1003)
assert candidate in metric and wheel in metric
noncand={k:v for k,v in metric.items() if k not in (candidate,wheel)}
assert len(noncand)==53
# Old trees must index exactly the same 53 topologies.
oldkeys=[(int(e['B']),int(e['I']),int(e['idx'])) for e in W.RALL['entries']]
assert len(oldkeys)==53 and set(oldkeys)==set(noncand)
mn=None;totleaf=totempty=totnode=totstress=0;t=time.time();reports=[]
for k in range(a.start,min(a.end,53)):
 entry=W.RALL['entries'][k];B,I,idx=oldkeys[k]
 fs=glob.glob(f'/mnt/data/resexact_{k:02d}_{B}_{idx}.json');assert len(fs)==1,(k,fs)
 C=json.load(open(fs[0]));assert int(C['entry_index'])==k and (int(C['B']),int(C['I']),int(C['idx']))==(B,I,idx)
 faces=[tuple(map(int,f)) for f in C['faces']];assert tuple(faces)==noncand[(B,I,idx)]
 assert int(C['qwb'])==W.QWBITS and int(C['qdb'])==W.QDBITS and int(C['lmixb'])==W.LMIXBITS
 pairs=[tuple(map(int,p)) for p in C['pairs']];rows=[]
 for sr in C['stresses']:
  nums=list(map(int,sr['weight_nums']));assert len(nums)==len(W.graph(B,I,faces)[0]) and min(nums)>=0 and sum(nums)==W.QW
  pp,dn,ex=W.exact_D(B,I,faces,nums);assert pp==pairs and dn==list(map(int,sr['Dnums']))
  # Explicit downward validity of every stored coefficient.
  for n,z in zip(dn,ex):assert F(n,W.QD)<=z
  rows.append([F(n,W.QD) for n in dn])
 oldleaves=W.reconstruct_old(entry);oldmap={int(o):(lo,hi) for o,lo,hi in oldleaves}
 roots=C['roots'];assert len(roots)==len(oldleaves) and {int(r['old_node']) for r in roots}==set(oldmap)
 lm=F(0) if False else None;nl=ne=0
 for r in roots:
  x,y,z=W.verify_tree(*oldmap[int(r['old_node'])],r['nodes'],pairs,rows)
  nl+=x;ne+=y;lm=z if lm is None or (z is not None and z<lm) else lm
 assert nl+ne>0
 mn=lm if mn is None or (lm is not None and lm<mn) else mn
 totnode+=sum(len(r['nodes']) for r in roots);totleaf+=nl;totempty+=ne;totstress+=len(rows)
 rec={'entry':k,'B':B,'I':I,'idx':idx,'stresses':len(rows),'old_leaves':len(oldleaves),'nodes':sum(len(r['nodes']) for r in roots),'leaves':nl,'empty':ne,'minmargin':float(lm)};reports.append(rec)
 print('VERIFIED',json.dumps(rec), 'elapsed',time.time()-t,flush=True)
out={'start':a.start,'end':min(a.end,53),'entries':len(reports),'stresses':totstress,'nodes':totnode,'leaves':totleaf,'empty':totempty,'minmargin':None if mn is None else [str(mn.numerator),str(mn.denominator)],'minmargin_float':None if mn is None else float(mn),'reports':reports,'sec':time.time()-t}
if a.out:open(a.out,'w').write(json.dumps(out,indent=2))
print('RESIDUAL CHUNK VERIFIED',json.dumps({k:out[k] for k in ['start','end','entries','stresses','nodes','leaves','empty','minmargin_float','sec']}))
