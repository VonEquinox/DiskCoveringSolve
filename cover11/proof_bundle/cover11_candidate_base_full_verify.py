#!/usr/bin/env python3
"""Independent exact replay of the 100834 direct candidate leaves.
Deferred leaves are counted and handled by cover11_candidate_refined_full_verify.py.
"""
from fractions import Fraction as F
import argparse,json,glob,time,sys
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int);ap.add_argument('--out');a=ap.parse_args()
pairs,rows=core.load_D_down()
base=[]
for f in sorted(glob.glob('/mnt/data/basewit_[0-4].json')):base+=json.load(open(f))['records']
base.sort(key=lambda r:int(r['rank']));assert len(base)==106136
cert=json.load(open('/mnt/data/cover11_candidate_cert_T14438.json'));nodes=cert['nodes']
stack=[(0,(F(0),)*9,(F(1),)*9)];boxmap={};leaf_order=[]
while stack:
 idx,lo,hi=stack.pop();n=nodes[idx]
 if 'split' in n:
  k=int(n['split']);m=(lo[k]+hi[k])/2;h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
  stack.append((int(n['child'][1]),tuple(l1),hi));stack.append((int(n['child'][0]),lo,tuple(h0)))
 elif 'leaf' in n:
  boxmap[idx]=(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi));leaf_order.append(idx)
 elif 'empty' in n:
  assert core.tighten_exact(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi)) is None
 else:raise AssertionError((idx,n))
assert [int(r['node']) for r in base]==leaf_order
end=len(base) if a.end is None else min(a.end,len(base));mn=None;npass=ndeferred=0;t=time.time()
for rank in range(a.start,end):
 r=base[rank];w=r['witness']
 if w is None or w.get('kind')=='fail':ndeferred+=1;continue
 lo,hi=boxmap[int(r['node'])];v=core.bound_exact(lo,hi,w,pairs,rows);assert v is not None
 mar=v-core.TARGET;assert mar>=0,(rank,r['node'],w.get('kind'),float(mar))
 mn=mar if mn is None or mar<mn else mn;npass+=1
 if (rank-a.start+1)%1000==0:
  print('ranks',rank+1,'pass',npass,'deferred',ndeferred,'minmargin',None if mn is None else float(mn),'sec',time.time()-t,flush=True)
  if (rank-a.start+1)%5000==0:core.line_exact.cache_clear()
out={'start':a.start,'end':end,'ranks':end-a.start,'pass':npass,'deferred':ndeferred,
     'minmargin':None if mn is None else [str(mn.numerator),str(mn.denominator)],
     'minmargin_float':None if mn is None else float(mn),'sec':time.time()-t}
if a.out:open(a.out,'w').write(json.dumps(out,indent=2))
print('CANDIDATE BASE CHUNK VERIFIED',json.dumps({k:out[k] for k in ['start','end','ranks','pass','deferred','minmargin_float','sec']}))
