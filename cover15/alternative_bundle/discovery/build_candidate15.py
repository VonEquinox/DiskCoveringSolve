"""Untrusted candidate-tree discovery; every force and local leaf checked exactly."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
import json,gzip,time,traceback
import force15 as G
import system15 as S
from candidate_global15 import local
R=Path(__file__).resolve().parent

def make():
 st=time.time();B=S.B;faces=S.FACES;M=G.bounds(B,faces);lo,hi=G.rootbox(M);solver=G.Solver(B,faces);nodes=[None];stack=[(0,lo,hi,0)];ndual=nlocal=0;depmax=0;mint=None;minl=None
 while stack:
  ni,lo,hi,dep=stack.pop();depmax=max(depmax,dep);box=G.tighten(M,lo,hi)
  if box is None:nodes[ni]={'kind':'EMPTY'};continue
  lo,hi=box;ml=local(lo,hi)
  if ml is not None and ml>0:
   nodes[ni]={'kind':'LOCAL'};nlocal+=1;minl=ml if minl is None else min(minl,ml);continue
  cert,x,status=solver.solve_slsqp(lo,hi,maxiter=100)
  if cert is None and status['radius']>float(G.RU)+1e-6 and not status['success']:cert,x,status=solver.solve_slsqp(lo,hi,maxiter=350)
  if cert is not None:
   cert['kind']='DUAL';nodes[ni]=cert;ndual+=1;md=G.check(B,solver.ed,lo,hi,cert);mint=md if mint is None else min(mint,md)
  else:
   if len(nodes)>200000 or dep>70:raise RuntimeError(('unresolved',ni,dep,status))
   widths=[h-l for l,h in zip(lo,hi)];axis=max(range(B-1),key=lambda j:widths[j]);mid=(lo[axis]+hi[axis])/2
   a=len(nodes);nodes.extend([None,None]);nodes[ni]={'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[a,a+1]}
   lhi=list(hi);lhi[axis]=mid;rlo=list(lo);rlo[axis]=mid;stack.extend([(a+1,rlo,list(hi),dep+1),(a,list(lo),lhi,dep+1)])
  if (ndual+nlocal)%100==0:print('progress','nodes',len(nodes),'stack',len(stack),'force',ndual,'local',nlocal,'depth',depmax,'sec',time.time()-st,flush=True)
 assert all(n is not None for n in nodes)
 out={'version':1,'faces':faces,'nodes':nodes,'summary':{'force':ndual,'local':nlocal,'maxdepth':depmax,'minimum_force_margin':str(mint),'minimum_local_margin':str(minl),'seconds':time.time()-st}}
 with gzip.open(R/'candidate15_tree.json.gz','wt',encoding='utf8') as f:json.dump(out,f,separators=(',',':'))
 print('CANDIDATE COMPLETE',len(nodes),out['summary'],flush=True)
if __name__=='__main__':make()
