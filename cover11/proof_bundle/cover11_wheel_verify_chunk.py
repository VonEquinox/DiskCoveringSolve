#!/usr/bin/env python3
from fractions import Fraction as F
import argparse,json,sys,time
sys.path.insert(0,'/mnt/data');import cover11_candidate_exact_core as core
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,required=True);ap.add_argument('--end',type=int,required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
D=json.load(open('/mnt/data/cover11_wheel_cert_T14438.json'));N=D['nodes'];B=10;DELTA=F(*D['delta']);TARGET=F(1443877286317993,10**16);A={int(k):F(*v) for k,v in D['acoef'].items()};pairs=[(i,j) for i in range(B) for j in range(i+1,B)];row=[2*A[min(j-i,B-(j-i))] for i,j in pairs]
R=json.load(open('/mnt/data/cover11_wheel_upgrade_refined.json'));rmap={int(z['old_idx']):z for z in R['roots']}
def key(lo,hi):return tuple(z for p in zip(lo,hi) for z in p)
def canonical(lo,hi):
 best=None;box=None;lo=tuple(lo);hi=tuple(hi)
 for rev in (0,1):
  ll=lo[::-1] if rev else lo;hh=hi[::-1] if rev else hi
  for k in range(B):
   l=ll[k:]+ll[:k];h=hh[k:]+hh[:k];q=key(l,h)
   if best is None or q<best:best=q;box=(l,h)
 return box
boxes=[None]*len(N);boxes[0]=canonical((F(0),)*B,(F(1),)*B);leaves=[]
for idx,n in enumerate(N):
 lo,hi=boxes[idx]
 if 'split'in n:
  k=int(n['split']);m=(lo[k]+hi[k])/2
  for side,j in enumerate(n['child']):
   l=list(lo);h=list(hi)
   if side==0:h[k]=m
   else:l[k]=m
   cb=canonical(tuple(l),tuple(h));j=int(j)
   if boxes[j] is None:boxes[j]=cb
   else:assert boxes[j]==cb
 elif 'leaf'in n:leaves.append(idx)
 elif 'empty'in n:assert core.tighten_exact([x*DELTA for x in lo],[x*DELTA for x in hi]) is None
end=min(a.end,len(leaves));mn=None;direct=refined=0;t=time.time()
for rank in range(a.start,end):
 idx=leaves[rank];lo0,hi0=boxes[idx];lo=tuple(x*DELTA for x in lo0);hi=tuple(x*DELTA for x in hi0);v=core.bound_exact(lo,hi,{'kind':'row','row':0},pairs,[row])
 if v is not None and v>=TARGET:
  direct+=1;mar=v-TARGET;mn=mar if mn is None or mar<mn else mn;continue
 assert idx in rmap,(rank,idx,float(v-TARGET) if v else None)
 z=rmap[idx];nodes=z['nodes'];bb=[None]*len(nodes);bb[0]=(lo,hi)
 for k,n in enumerate(nodes):
  l,h=bb[k]
  if 'split'in n:
   j=int(n['split']);m=F(int(n['mid'][0]),int(n['mid'][1]));h0=list(h);h0[j]=m;l1=list(l);l1[j]=m;bb[int(n['child'][0])]=(l,tuple(h0));bb[int(n['child'][1])]=(tuple(l1),h)
  elif 'leaf'in n:
   vv=core.bound_exact(l,h,{'kind':'row','row':0},pairs,[row]);assert vv is not None and vv>=TARGET
   mar=vv-TARGET;mn=mar if mn is None or mar<mn else mn;refined+=1
  elif 'empty'in n:assert core.tighten_exact(l,h) is None
  else:raise AssertionError((idx,k,n))
out={'start':a.start,'end':end,'total_leaves':len(leaves),'direct':direct,'refined_leaves':refined,'minmargin':None if mn is None else [mn.numerator,mn.denominator],'minmargin_float':None if mn is None else float(mn),'sec':time.time()-t}
open(a.out,'w').write(json.dumps(out,indent=2));print(json.dumps(out))
