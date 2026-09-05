"""Untrusted numerical proposal generator; every stored leaf is exact-checked."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
import json,gzip,time,traceback,sys,math
import numpy as np
from concurrent.futures import ProcessPoolExecutor,as_completed
import force_fast20 as G
import system20 as S
from candidate_global20 import local
P=Path(__file__).resolve().parent;OUT=P/'forest_chunks20';OUT.mkdir(exist_ok=True)
MAP=json.load(open(P/'candidate20_map.json'));MET=json.load(open(P/'cases20_meta.json'));CAND=(MAP['B'],MAP['idx']);CACHE={}
def masks(B):
 if B not in CACHE:CACHE[B]=np.load(P/f'masks20_B{B}.npy',mmap_mode='r')
 return CACHE[B]
def make_tree(B,faces,is_candidate=False):
 M=G.bounds(B,faces);lo,hi=G.rootbox(M);solver=G.Solver(B,faces);nodes=[None];stack=[(0,lo,hi,0)];depthmax=0;nd=nl=0;st=time.time()
 while stack:
  ni,lo,hi,dep=stack.pop();depthmax=max(dep,depthmax);box=G.tighten(M,lo,hi)
  if box is None:nodes[ni]={'kind':'EMPTY'};continue
  lo,hi=box
  if is_candidate:
   ml=local(tuple(lo),tuple(hi))
   if ml is not None and ml>0:nodes[ni]={'kind':'LOCAL'};nl+=1;continue
  cert,x,status=solver.solve_slsqp(lo,hi,maxiter=80)
  if cert is None and status['radius']>float(G.RU)+1e-6 and not status['success']:cert,x,status=solver.solve_slsqp(lo,hi,maxiter=350)
  if cert is not None:
   assert G.check(B,solver.ed,lo,hi,cert)>0;nodes[ni]={'kind':'DUAL','cert':cert};nd+=1
  else:
   if len(nodes)>50000 or dep>70:raise RuntimeError(('unresolved',ni,dep,status,lo,hi))
   widths=[h-l for l,h in zip(lo,hi)];q=np.array(x[:2*(B-1)]).reshape(-1,2);deficits=1-(q*q).sum(axis=1)
   axis=max(range(B-1),key=lambda j:max(0,float(deficits[j])) if widths[j]>F(1,10**10) else -1)
   if deficits[axis]<1e-10:axis=max(range(B-1),key=lambda j:widths[j])
   mid=(lo[axis]+hi[axis])/2;a=len(nodes);nodes.extend([None,None]);nodes[ni]={'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[a,a+1]}
   lhi=list(hi);lhi[axis]=mid;rlo=list(lo);rlo[axis]=mid
   stack.extend([(a+1,rlo,list(hi),dep+1),(a,list(lo),lhi,dep+1)])
  if is_candidate and (nd+nl)%200==0:print('candidate progress',len(nodes),'depth',depthmax,'dual',nd,'local',nl,'sec',time.time()-st,flush=True)
 assert all(n is not None for n in nodes)
 return nodes,depthmax

def chunk(task):
 B,lo,hi=task;path=OUT/f'B{B}_{lo:06}_{hi:06}.jsonl.gz';tmp=Path(str(path)+'.tmp');st=time.time();total=0;failed=[]
 if path.exists():return task,'cached',0,0
 try:
  with gzip.open(tmp,'wt',encoding='utf8',compresslevel=6) as f:
   for j in range(lo,hi):
    if (B,j)==CAND:continue
    fa=[[v for v in range(20) if int(m)>>v&1] for m in masks(B)[j]]
    try:
     nodes,dep=make_tree(B,fa)
     f.write(json.dumps({'B':B,'idx':j,'nodes':nodes},separators=(',',':'))+'\n');total+=len(nodes)
    except Exception as e:
     failed.append({'B':B,'idx':j,'error':repr(e),'traceback':traceback.format_exc()})
  if failed:
   (OUT/f'FAILED_B{B}_{lo}.json').write_text(json.dumps(failed,indent=2));tmp.unlink(missing_ok=True);return task,'failed',total,time.time()-st
  tmp.rename(path);return task,'done',total,time.time()-st
 except Exception as e:
  return task,'fatal '+repr(e),total,time.time()-st

if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='candidate':
  st=time.time();nodes,dep=make_tree(S.B,S.FACES,True)
  with gzip.open(P/'candidate20_tree.json.gz','wt',encoding='utf8') as f:json.dump({'version':1,'faces':S.FACES,'nodes':nodes},f,separators=(',',':'))
  print('CANDIDATE COMPLETE',len(nodes),'maxdepth',dep,'seconds',time.time()-st,flush=True)
 else:
  workers=int(sys.argv[2]) if len(sys.argv)>2 else 4;tasks=[];step=500
  for family in MET['families']:
   B=family['B'];n=family['orbits']
   for i in range(0,n,step):tasks.append((B,i,min(n,i+step)))
  if mode=='test':tasks=[(12,0,50),(13,0,50),(14,0,50),(15,0,50),(16,0,50)]
  st=time.time();count=cases=nn=0;fails=[]
  with ProcessPoolExecutor(workers) as ex:
   fs=[ex.submit(chunk,t) for t in tasks]
   for fu in as_completed(fs):
    task,status,n,t=fu.result();count+=1;cases+=task[2]-task[1];nn+=n
    if status.startswith(('failed','fatal')):fails.append(task)
    if count%10==0 or status.startswith(('failed','fatal')) or mode=='test':print('CHUNKS',count,'/',len(tasks),'cases',cases,'nodes',nn,'failed',fails,'sec',time.time()-st,'last',task,status,t,flush=True)
  print('FOREST DONE',count,'cases',cases,'nodes',nn,'failed',fails,'seconds',time.time()-st,flush=True)
