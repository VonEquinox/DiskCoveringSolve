#!/usr/bin/env python3
"""Exact integer-flow replay for every noncandidate residual triangulation.
All acceptance uses standard-library exact arithmetic.
"""
from pathlib import Path
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
from fractions import Fraction as F
from math import isqrt
from collections import Counter
import json,gzip,hashlib,time,sys
from exact_arcs import QA,QB,arc_halfplane
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'enumeration'
sys.path.insert(0,str(ROOT.parent))
from verify_partition import signature,A,C,D
QL=10**15;DF=10**36;RU=F(346645456927389644,10**18)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilroot(n):
 assert n>=0;r=isqrt(n);return r+(r*r<n)
def edges(B,faces):
 assert len(faces)==24-B
 return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+13+k,B+c) for k,fa in enumerate(faces) for c in fa]
def bounds(B,faces):
 M=[[10**8]*B for _ in range(B)]
 for i in range(B):M[i][i]=0
 for i in range(B-1):M[i][i+1]=A;M[i+1][i]=0
 M[B-1][0]=A-D;M[0][B-1]=D
 key=signature(B,faces);assert key is not None
 for mask in key:
  start=[i for i in range(B) if mask>>i&1 and not(mask>>((i-1)%B)&1)]
  end=[i for i in range(B) if not(mask>>i&1) and mask>>((i-1)%B)&1]
  assert len(start)==len(end)==1;i,j=start[0],end[0]
  M[i][j]=min(M[i][j],C-(D if j<i else 0))
 for k in range(B):
  for i in range(B):
   for j in range(B):M[i][j]=min(M[i][j],M[i][k]+M[k][j])
 assert all(M[i][i]>=0 for i in range(B))
 return M

def rootbox(M):return [F(-M[i][0],D) for i in range(1,len(M))],[F(M[0][i],D) for i in range(1,len(M))]
def tighten(M,lo,hi):
 B=len(M);l=[F(0)]+list(lo);h=[F(0)]+list(hi)
 nh=[min(h[i]+F(M[i][j],D) for i in range(B)) for j in range(B)]
 nl=[max(l[i]-F(M[j][i],D) for i in range(B)) for j in range(B)]
 if any(a>b for a,b in zip(nl,nh)):return None
 assert nl[0]<=0<=nh[0]
 return nl[1:],nh[1:]

def check_dual(B,ed,lo,hi,r):
 lam=r['lambda_num'];eta=r['eta_num'];dv=r['support_vector_num'];flow=r['force_num']
 assert len(lam)==len(eta)==len(dv)==B-1
 assert all(type(x)is int and x>=0 for x in lam+eta)
 assert all(len(z)==2 and all(type(v)is int for v in z) for z in dv+flow)
 assert len(flow)==len(ed)
 aa=[];bb=[]
 for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
 target=[[0,0] for _ in range(37)]
 for i in range(B-1):
  for d in range(2):
   target[i+1][d]=(lam[i]*aa[i][d]-eta[i]*dv[i][d])*(DF//(QL*QA))
   target[0][d]-=target[i+1][d]
 act=[[0,0] for _ in range(37)]
 for (u,v),z in zip(ed,flow):
  for d in range(2):act[u][d]+=z[d];act[v][d]-=z[d]
 assert act==target,'integer equilibrium failed'
 ns=sum(ceilroot(z[0]**2+z[1]**2) for z in flow)
 bn=sum(lam[i]*(bb[i]-aa[i][0]*(QB//QA))-eta[i]*(ceilroot(dv[i][0]**2+dv[i][1]**2)-dv[i][0])*(QB//QA) for i in range(B-1))
 assert ns>0 and bn>0
 assert bn*RU.denominator*DF>RU.numerator*ns*QL*QB,('dual failed',float(F(bn*DF,QL*QB*ns)-RU))
 return F(bn*DF,QL*QB*ns)-RU

def verify():
 st=time.time();rr=json.loads((OUT/'metric_residuals_new.json').read_text());cm=json.loads((OUT/'candidate_orbit_map.json').read_text());cand=cm['residual_index']
 cp=OUT/'noncandidate_certificate.json.gz'
 with gzip.open(cp,'rt',encoding='utf8') as f:cert=json.load(f)
 assert cert['schema']=='all-noncandidate-integer-force-v1' and F(cert['radius_upper'])==RU and cert['lambda_scale']==QL and cert['flow_scale']==DF
 assert cert['residual_sha256']==sha(OUT/'metric_residuals_new.json')
 covered=set();counts=Counter();tot=0;mindual=None;maxdepth=0
 for case in cert['cases']:
  k=case['residual_index'];assert type(k)is int and 0<=k<len(rr) and k!=cand and k not in covered;covered.add(k)
  row=rr[k]
  for key in ('B','I','idx','faces'):assert case[key]==row[key]
  B=row['B'];ed=edges(B,row['faces']);M=bounds(B,row['faces']);lo,hi=rootbox(M)
  ns=case['nodes'];seen=set();stack=[(0,lo,hi,0)];tot+=len(ns)
  while stack:
   idx,lo,hi,depth=stack.pop();assert type(idx)is int and 0<=idx<len(ns) and idx not in seen;seen.add(idx);maxdepth=max(maxdepth,depth)
   node=ns[idx];kind=node['kind'];box=tighten(M,lo,hi)
   if kind=='EMPTY':assert box is None;counts[kind]+=1;continue
   assert box is not None;lo,hi=box
   if kind=='SPLIT':
    j=node['axis'];m=F(node['mid']);children=node['children']
    assert type(j)is int and 0<=j<B-1 and lo[j]<m<hi[j] and len(children)==2 and children[0]!=children[1]
    lh=list(hi);lh[j]=m;rl=list(lo);rl[j]=m
    stack.extend([(children[1],rl,list(hi),depth+1),(children[0],list(lo),lh,depth+1)])
   elif kind=='DUAL':
    mar=check_dual(B,ed,lo,hi,node);mindual=mar if mindual is None else min(mindual,mar);counts[kind]+=1
   else:raise AssertionError(('bad node',kind))
  assert seen==set(range(len(ns)))
 assert covered==set(range(len(rr)))-{cand}
 out={'verified':True,'cases':len(covered),'total_nodes':tot,'leaf_counts':dict(counts),'maxdepth':maxdepth,'minimum_radius_margin':str(mindual),'minimum_radius_margin_float':float(mindual),'certificate_sha256':sha(cp),'seconds':time.time()-st}
 (OUT/'noncandidate_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 verify()
