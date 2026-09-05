"""Numerical discovery; each proposed leaf must pass exact local/force checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
import json,gzip,time,traceback
import force15 as G
import candidate_global15 as L
import system15 as S
R=Path(__file__).resolve().parent

def make():
 st=time.time();B=S.B;faces=[list(f) for f in S.FACES];M=G.bounds(B,faces);lo,hi=G.rootbox(M);solver=G.Solver(B,faces);nodes=[None];stack=[(0,lo,hi,0)];nd=nl=ne=0;minmar=None;minlocal=None;maxdepth=0
 while stack:
  ni,lo,hi,dep=stack.pop();box=G.tighten(M,lo,hi);maxdepth=max(maxdepth,dep)
  if box is None:nodes[ni]={'kind':'EMPTY'};ne+=1;continue
  lo,hi=box;l=L.local(lo,hi)
  if l is not None and l>=0:
   nodes[ni]={'kind':'LOCAL'};nl+=1;minlocal=l if minlocal is None else min(minlocal,l);continue
  cert,x,status=solver.solve_slsqp(lo,hi,maxiter=100)
  if cert is None and status['radius']>float(G.RU)+1e-6 and (not status['success'] or status.get('exact_margin',-1)<0):cert,x,status=solver.solve_slsqp(lo,hi,maxiter=450)
  if cert is not None:
   cert['kind']='DUAL';nodes[ni]=cert;nd+=1;mar=G.check(B,solver.ed,lo,hi,cert);minmar=mar if minmar is None else min(minmar,mar);continue
  if dep>70 or len(nodes)>100000:raise RuntimeError(('unresolved',ni,dep,status,[str(x) for x in lo],[str(x) for x in hi]))
  widths=[h-l for l,h in zip(lo,hi)];axis=max(range(B-1),key=lambda i:widths[i]);mid=(lo[axis]+hi[axis])/2
  assert lo[axis]<mid<hi[axis]
  a=len(nodes);nodes.extend([None,None]);nodes[ni]={'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[a,a+1]}
  lhi=list(hi);lhi[axis]=mid;rlo=list(lo);rlo[axis]=mid
  stack.extend([(a+1,rlo,list(hi),dep+1),(a,list(lo),lhi,dep+1)])
  if len(nodes)%100==1:print('nodes',len(nodes),'dual',nd,'local',nl,'empty',ne,'depth',maxdepth,'sec',time.time()-st,flush=True)
 out={'version':1,'B':B,'I':S.N-B,'faces':faces,'bounds_matrix':M,'nodes':nodes,'summary':{'nodes':len(nodes),'dual':nd,'local':nl,'empty':ne,'maxdepth':maxdepth,'min_force_margin':str(minmar),'min_local_margin':str(minlocal),'seconds':time.time()-st}}
 with gzip.open(R/'candidate15_tree.json.gz','wt',encoding='utf8') as f:json.dump(out,f,separators=(',',':'))
 print('COMPLETE',json.dumps(out['summary']),flush=True)
if __name__=='__main__':
 try:make()
 except Exception as e:
  (R/'candidate15_tree.failed.txt').write_text(traceback.format_exc());raise
