#!/usr/bin/env python3
"""Discover exact force trees for all noncandidate residuals. Explicitly fail on pending leaves."""
from pathlib import Path
import json,gzip,time,os,traceback
from multiprocessing import Pool
from fractions import Fraction as F
import numpy as np
from force_solver import Solver,QL,DF,RU
from angle_domains import bounds,rootbox,tighten
BASE=Path(__file__).resolve().parent;ENUM=BASE/'enumeration';CASES=BASE/'noncandidate_cases';CASES.mkdir(exist_ok=True)
ROWS=json.loads((ENUM/'metric_residuals.json').read_text());CAND=685

def one(k):
 dst=CASES/f'{k}.json.gz'
 if dst.exists():
  with gzip.open(dst,'rt') as f:c=json.load(f)
  if c.get('complete'):return k,len(c['nodes']),c.get('min_margin',0),True
 st=time.time();row=ROWS[k];B=row['B'];faces=row['faces'];s=Solver(B,faces);M=bounds(B,faces);lo,hi=rootbox(M)
 nodes=[None];stack=[(0,lo,hi,None,0)];nsolve=0;mindual=1.;maxdepth=0
 try:
  while stack:
   idx,lo,hi,x0,dep=stack.pop();box=tighten(M,lo,hi);maxdepth=max(maxdepth,dep)
   if box is None:nodes[idx]={'kind':'EMPTY'};continue
   lo,hi=box
   mar,rec,opt,z=s.solve(lo,hi,x0);nsolve+=1
   if mar>0:nodes[idx]=rec;mindual=min(mindual,float(mar));continue
   if nsolve>20000 or dep>70:
    nodes[idx]={'kind':'PENDING'}
    raise RuntimeError(('budget',nsolve,dep,float(mar),opt.fun))
   widths=np.array([float(h-l) for l,h in zip(lo,hi)])
   lm=np.array(rec['lambda_num'],float)/QL
   scores=(lm+max(lm.max(),.001)*.05)*widths**2
   j=int(np.argmax(scores));mid=(lo[j]+hi[j])/2;assert lo[j]<mid<hi[j]
   ch=[len(nodes),len(nodes)+1];nodes.extend([None,None]);nodes[idx]={'kind':'SPLIT','axis':j,'mid':str(mid),'children':ch}
   lh=list(hi);lh[j]=mid;rl=list(lo);rl[j]=mid
   stack.extend([(ch[1],rl,list(hi),opt.x,dep+1),(ch[0],list(lo),lh,opt.x,dep+1)])
  out={'complete':True,'residual_index':k,**{a:row[a] for a in ('B','I','idx','faces')},'nodes':nodes,'min_margin':mindual,'maxdepth':maxdepth,'seconds':time.time()-st}
  with gzip.open(dst,'wt') as f:json.dump(out,f,separators=(',',':'))
  return k,len(nodes),mindual,True
 except Exception as e:
  out={'complete':False,'residual_index':k,'error':repr(e),'trace':traceback.format_exc(),'nodes':nodes,'remaining_stack':len(stack),'seconds':time.time()-st}
  with gzip.open(dst,'wt') as f:json.dump(out,f,separators=(',',':'))
  return k,len(nodes),0.,False

def main():
 st=time.time();done=0;total=0;bad=[];small=1.;hard=[]
 with Pool(4) as pool:
  for k,nn,mar,ok in pool.imap_unordered(one,[k for k in range(len(ROWS)) if k!=CAND],chunksize=1):
   done+=1;total+=nn
   if not ok:bad.append(k)
   else:small=min(small,mar)
   if nn>20:hard.append([k,nn,ok]);print('hard',k,nn,ok,flush=True)
   if done%200==0 or not ok:
    report={'done':done,'total_nodes':total,'failures':bad,'hard':hard,'min_radius_margin_float':small,'seconds':time.time()-st};(BASE/'noncandidate_progress.json').write_text(json.dumps(report,indent=2));print(done,total,len(bad),round(time.time()-st,1),flush=True)
 out={'done':done,'total_nodes':total,'failures':bad,'hard':hard,'min_radius_margin_float':small,'seconds':time.time()-st,'complete':done==len(ROWS)-1 and not bad}
 (BASE/'noncandidate_progress.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
