#!/usr/bin/env python3
"""Generate a full candidate-angle partition, DUAL or rigorously LOCAL at every leaf."""
from pathlib import Path
import json,gzip,time,traceback
from fractions import Fraction as F
import numpy as np
from framework import ALL_FACES
from force_solver import Solver,QL
from angle_domains import bounds,rootbox,tighten
from local_box import accepts
BASE=Path(__file__).resolve().parent

def main():
 st=time.time();B=10;s=Solver(B,ALL_FACES);M=bounds(B,ALL_FACES);lo,hi=rootbox(M)
 nodes=[None];stack=[(0,lo,hi,None,0)];nsolve=0;mindual=1.;maxdepth=0;counts={'DUAL':0,'LOCAL':0,'EMPTY':0,'SPLIT':0}
 try:
  while stack:
   idx,lo,hi,x0,dep=stack.pop();box=tighten(M,lo,hi);maxdepth=max(maxdepth,dep)
   if box is None:nodes[idx]={'kind':'EMPTY'};counts['EMPTY']+=1;continue
   lo,hi=box
   if accepts(lo,hi):nodes[idx]={'kind':'LOCAL'};counts['LOCAL']+=1;continue
   mar,rec,opt,z=s.solve(lo,hi,x0);nsolve+=1
   if mar>0:nodes[idx]=rec;counts['DUAL']+=1;mindual=min(mindual,float(mar));continue
   if nsolve>100000 or dep>100:raise RuntimeError(('budget',nsolve,dep,float(mar),opt.fun))
   widths=np.array([float(h-l) for l,h in zip(lo,hi)]);lm=np.array(rec['lambda_num'],float)/QL
   scores=(lm+max(lm.max(),.001)*.05)*widths**2;j=int(np.argmax(scores));mid=(lo[j]+hi[j])/2;assert lo[j]<mid<hi[j]
   ch=[len(nodes),len(nodes)+1];nodes.extend([None,None]);nodes[idx]={'kind':'SPLIT','axis':j,'mid':str(mid),'children':ch};counts['SPLIT']+=1
   lh=list(hi);lh[j]=mid;rl=list(lo);rl[j]=mid
   stack.extend([(ch[1],rl,list(hi),opt.x,dep+1),(ch[0],list(lo),lh,opt.x,dep+1)])
   if nsolve%200==0:
    print(nsolve,len(nodes),len(stack),counts,round(time.time()-st,1),flush=True)
  out={'complete':True,'B':10,'I':4,'faces':ALL_FACES,'residual_index':685,'nodes':nodes,'counts':counts,'min_margin':mindual,'maxdepth':maxdepth,'seconds':time.time()-st}
 except Exception as e:out={'complete':False,'nodes':nodes,'error':repr(e),'trace':traceback.format_exc(),'pending':len(stack),'counts':counts}
 with gzip.open(BASE/'candidate_certificate.json.gz','wt') as f:json.dump(out,f,separators=(',',':'))
 print(json.dumps({k:v for k,v in out.items() if k!='nodes'},indent=2),flush=True)
if __name__=='__main__':main()
