"""Discovery only. Every finished case will be independently replayed."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
import json,gzip,time,traceback,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
import force15 as G
R=Path(__file__).resolve().parent;FAMILY=int(sys.argv[1]);OUT=R/f'forest15_B{FAMILY}';OUT.mkdir(exist_ok=True)
ROWS=json.load(open(R/f'metric15_residuals_B{FAMILY}.json'));CAND=json.load(open(R/'candidate15_map_family.json'))['family_residual_index'] if FAMILY==11 else -1

def run_case(k):
 path=OUT/f'{k:05}.json.gz'
 if path.exists():return k,'cached',0,0
 st=time.time();row=ROWS[k];B=row['B'];faces=row['faces'];M=G.bounds(B,faces);lo,hi=G.rootbox(M);solver=G.Solver(B,faces);nodes=[None];stack=[(0,lo,hi,0)];marmin=None;depthmax=0
 try:
  while stack:
   ni,lo,hi,dep=stack.pop();box=G.tighten(M,lo,hi);depthmax=max(depthmax,dep)
   if box is None:nodes[ni]={'kind':'EMPTY'};continue
   lo,hi=box
   cert,x,status=solver.solve_slsqp(lo,hi,maxiter=70)
   if cert is None and status['radius']>float(G.RU)+1e-5 and (not status['success'] or status.get('exact_margin',-1)<0):cert,x,status=solver.solve_slsqp(lo,hi,maxiter=300)
   if cert is not None:
    cert['kind']='DUAL';nodes[ni]=cert;mar=G.check(B,solver.ed,lo,hi,cert);marmin=mar if marmin is None else min(marmin,mar);continue
   if len(nodes)>10000 or dep>50:raise RuntimeError(('unresolved',ni,dep,status,[float(v) for v in lo],[float(v) for v in hi]))
   widths=[h-l for l,h in zip(lo,hi)];axis=max(range(B-1),key=lambda j:widths[j]);mid=(lo[axis]+hi[axis])/2
   a=len(nodes);nodes.extend([None,None]);nodes[ni]={'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[a,a+1]}
   lhi=list(hi);lhi[axis]=mid;rlo=list(lo);rlo[axis]=mid
   stack.extend([(a+1,rlo,list(hi),dep+1),(a,list(lo),lhi,dep+1)])
  out={'version':1,'residual_index':k,**row,'nodes':nodes,'summary':{'nodes':len(nodes),'maxdepth':depthmax,'min_margin':str(marmin),'seconds':time.time()-st}}
  with gzip.open(path,'wt',encoding='utf8') as f:json.dump(out,f,separators=(',',':'))
  return k,'done',len(nodes),time.time()-st
 except Exception as e:
  data={'case':k,'error':repr(e),'traceback':traceback.format_exc(),'nodes':len(nodes),'seconds':time.time()-st};(OUT/f'{k:05}.failed.json').write_text(json.dumps(data,indent=2));return k,'failed',len(nodes),time.time()-st

if __name__=='__main__':
 workers=int(sys.argv[2]) if len(sys.argv)>2 else 4;st=time.time();count=0;totalnodes=0;failed=[]
 with ProcessPoolExecutor(max_workers=workers) as ex:
  tasks=[ex.submit(run_case,k) for k in range(len(ROWS)) if k!=CAND]
  for fut in as_completed(tasks):
   k,status,n,sec=fut.result();count+=1;totalnodes+=n
   if status=='failed':failed.append(k);print('FAILED',k,'nodes',n,'sec',sec,flush=True)
   if count%100==0 or status=='failed':print('progress',count,'/',len(tasks),'nodes',totalnodes,'failed',len(failed),'seconds',time.time()-st,flush=True)
 print('DONE',count,'failed',failed,'seconds',time.time()-st,flush=True)
