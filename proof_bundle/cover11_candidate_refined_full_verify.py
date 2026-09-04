#!/usr/bin/env python3
"""Independent exact replay of the 5302 refined candidate roots.

External leaves are checked by exact graph-energy affine minorants. Local
leaves are checked against the exact full-KKT gap center using the rigorous
1e-29 gap-coordinate transfer bound.
"""
from fractions import Fraction as F
import argparse,glob,json,sys,time
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int);ap.add_argument('--out');a=ap.parse_args()
pairs,rows=core.load_D_down()
# exact old leaf boxes
old=json.load(open('/mnt/data/cover11_candidate_cert_T14438.json'))['nodes']
stack=[(0,(F(0),)*9,(F(1),)*9)];oldbox={};leaf_order=[]
while stack:
 idx,lo,hi=stack.pop();n=old[idx]
 if 'split' in n:
  k=int(n['split']);m=(lo[k]+hi[k])/2;h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
  stack.append((int(n['child'][1]),tuple(l1),hi));stack.append((int(n['child'][0]),lo,tuple(h0)))
 elif 'leaf' in n:
  oldbox[idx]=(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi));leaf_order.append(idx)
 elif 'empty' in n:assert core.tighten_exact(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi)) is None
# identify exact deferred set from base records
base=[]
for f in sorted(glob.glob('/mnt/data/basewit_[0-4].json')):base+=json.load(open(f))['records']
base.sort(key=lambda r:int(r['rank']));assert [int(r['rank']) for r in base]==list(range(106136));assert [int(r['node']) for r in base]==leaf_order
deferred=sorted(int(r['node']) for r in base if r['witness'] is None or r['witness'].get('kind')=='fail')
assert len(deferred)==5302
roots=json.load(open('/mnt/data/cover11_candidate_refined_rational_splits.json'))['roots'];rmap={int(r['old_node']):r for r in roots};assert set(deferred)<=set(rmap)
# candidate center and proven transfer error
L=json.load(open('/mnt/data/cover11_local_regularized_certificate.json'));center=[F(x) for x in L['gap_center']]
# ninth center gap comes from the high-precision exact-root proposal and is
# independently within EPS by cover11_local_theorem_verify_v2.py.
H=json.load(open('/mnt/data/cover11_highprec_stresses.json'));center.append(F(H['gaps'][8]))
EPS=F(1,10**29);RLOC=F(11,5000)
end=len(deferred) if a.end is None else min(a.end,len(deferred));mn=None;nexternal=nlocal=nempty=nnodes=0;t=time.time()
for rank in range(a.start,end):
 oldnode=deferred[rank];ns=rmap[oldnode]['nodes'];boxes=[None]*len(ns);boxes[0]=oldbox[oldnode];nnodes+=len(ns)
 for idx,n in enumerate(ns):
  assert boxes[idx] is not None,(rank,oldnode,idx)
  lo,hi=boxes[idx]
  if 'split' in n:
   k=int(n['split']);m=F(int(n['mid'][0]),int(n['mid'][1]));assert lo[k]<m<hi[k]
   h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m;c0,c1=map(int,n['child']);assert boxes[c0] is None and boxes[c1] is None
   boxes[c0]=(lo,tuple(h0));boxes[c1]=(tuple(l1),hi)
  elif n.get('kind') in ('row','mix'):
   v=core.bound_exact(lo,hi,n,pairs,rows);assert v is not None;mar=v-core.TARGET;assert mar>=0,(rank,oldnode,idx,float(mar))
   mn=mar if mn is None or mar<mn else mn;nexternal+=1
  elif 'local' in n:
   th=core.tighten_exact(lo,hi);assert th is not None;l,h=th;d2=F(0)
   for j in range(9):
    far=max(abs(l[j]-(center[j]-EPS)),abs(l[j]-(center[j]+EPS)),abs(h[j]-(center[j]-EPS)),abs(h[j]-(center[j]+EPS)))
    d2+=far*far
   assert d2<RLOC*RLOC,(rank,oldnode,idx,float(d2-RLOC*RLOC));nlocal+=1
  elif 'empty' in n:
   assert core.tighten_exact(lo,hi) is None;nempty+=1
  else:raise AssertionError((rank,oldnode,idx,n))
 if (rank-a.start+1)%100==0:
  print('roots',rank+1,'external',nexternal,'local',nlocal,'minmargin',None if mn is None else float(mn),'sec',time.time()-t,flush=True)
  core.line_exact.cache_clear()
out={'start':a.start,'end':end,'roots':end-a.start,'nodes':nnodes,'external':nexternal,'local':nlocal,'empty':nempty,'minmargin':None if mn is None else [mn.numerator,mn.denominator],'minmargin_float':None if mn is None else float(mn),'sec':time.time()-t}
if a.out:open(a.out,'w').write(json.dumps(out,indent=2))
print('CANDIDATE REFINED CHUNK VERIFIED',json.dumps({k:out[k] for k in ['start','end','roots','nodes','external','local','empty','minmargin_float','sec']}))
