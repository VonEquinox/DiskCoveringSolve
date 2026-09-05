#!/usr/bin/env python3
"""Exact integer-force inequalities and exhaustive tree replay. No numerical optimizer imports."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from collections import Counter
import json,gzip,time,hashlib
from exact_arcs import QA,QB,arc_halfplane
from angle_domains import bounds,rootbox,tighten
from verify_partition import RU
from framework import ALL_FACES,ACTIVE_FACES
from local_box import bound as local_bound
from anchor_isolation import RADIUS
if not __debug__:raise RuntimeError('Assertions must be enabled')
BASE=Path(__file__).resolve().parent;ENUM=BASE/'enumeration'
QL=10**15;DF=10**36;NN=40

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def ceilroot(n):
 assert type(n)is int and n>=0;r=isqrt(n);return r+(r*r<n)
def edges(B,faces):
 assert len(faces)==26-B
 return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+14+k,B+c) for k,fa in enumerate(faces) for c in fa]

def check_dual(B,ed,lo,hi,r):
 lam=r['lambda_num'];eta=r['eta_num'];dv=r['support_vector_num'];flow=r['force_num']
 assert len(lam)==len(eta)==len(dv)==B-1 and len(flow)==len(ed)
 assert all(type(x)is int and x>=0 for x in lam+eta)
 assert all(len(z)==2 and all(type(v)is int for v in z) for z in dv+flow)
 assert DF%(QL*QA)==0 and QB%QA==0
 planes=[arc_halfplane(l,h) for l,h in zip(lo,hi)];target=[[0,0] for _ in range(NN)]
 for i,((a,b),dd) in enumerate(zip(planes,dv)):
  for d in range(2):
   target[i+1][d]=(lam[i]*a[d]-eta[i]*dd[d])*(DF//(QL*QA));target[0][d]-=target[i+1][d]
 actual=[[0,0] for _ in range(NN)]
 for (u,v),f in zip(ed,flow):
  assert 0<=u<NN and 0<=v<NN and u!=v
  for d in range(2):actual[u][d]+=f[d];actual[v][d]-=f[d]
 assert actual==target,'exact integer equilibrium failed'
 ns=sum(ceilroot(f[0]**2+f[1]**2) for f in flow)
 bn=sum(lam[i]*(b-a[0]*(QB//QA))-eta[i]*(ceilroot(dd[0]**2+dd[1]**2)-dd[0])*(QB//QA) for i,((a,b),dd) in enumerate(zip(planes,dv)))
 assert ns>0 and bn>0
 assert bn*DF*RU.denominator>RU.numerator*ns*QL*QB,'strict force inequality failed'
 return F(bn*DF,QL*QB*ns)-RU

def replay_tree(case,allow_local=False):
 B=case['B'];faces=case['faces'];ed=edges(B,faces);M=bounds(B,faces);lo,hi=rootbox(M)
 if allow_local:
  assert B==10 and faces==[list(f) for f in ALL_FACES]
  # Selecting these face-witness nodes embeds the 44-rod local graph in the full graph.
  full={tuple(sorted(e)) for e in ed}
  for fa in ACTIVE_FACES:
   k=faces.index(list(fa))
   for c in fa:assert tuple(sorted((B+14+k,B+c))) in full
 ns=case['nodes'];assert isinstance(ns,list) and ns
 seen=set();stack=[(0,lo,hi,0)];counts=Counter();mindual=None;maxlocal=F(0);maxdepth=0
 while stack:
  i,lo,hi,dep=stack.pop();assert type(i)is int and 0<=i<len(ns) and i not in seen;seen.add(i);maxdepth=max(maxdepth,dep)
  r=ns[i];assert isinstance(r,dict);kind=r['kind'];box=tighten(M,lo,hi)
  if kind=='EMPTY':assert box is None;counts[kind]+=1;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   j=r['axis'];mid=F(r['mid']);ch=r['children'];assert type(j)is int and 0<=j<B-1 and lo[j]<mid<hi[j]
   assert isinstance(ch,list) and len(ch)==2 and ch[0]!=ch[1]
   lh=list(hi);lh[j]=mid;rl=list(lo);rl[j]=mid
   stack.extend([(ch[1],rl,list(hi),dep+1),(ch[0],list(lo),lh,dep+1)])
  elif kind=='DUAL':
   mar=check_dual(B,ed,lo,hi,r);mindual=mar if mindual is None else min(mindual,mar)
  elif kind=='LOCAL':
   assert allow_local;v=local_bound(lo,hi);assert v<=RADIUS**2;maxlocal=max(maxlocal,v)
  else:raise AssertionError(('unproved leaf',kind))
  counts[kind]+=1
 assert seen==set(range(len(ns)))
 assert counts['SPLIT']*2+1==len(ns)
 return {'counts':counts,'nodes':len(ns),'maxdepth':maxdepth,'mindual':mindual,'maxlocal':maxlocal}

def check_residual_partition(keys,n):
 assert all(type(k)is int for k in keys)
 assert len(keys)==len(set(keys)) and set(keys)==set(range(n))-{685},'incomplete residual partition'

def verify_noncandidate():
 st=time.time();rr=json.loads((ENUM/'metric_residuals.json').read_text());path=BASE/'noncandidate_certificate.json.gz'
 with gzip.open(path,'rt') as f:cert=json.load(f)
 assert cert['schema']=='r14-all-noncandidate-force-v1' and F(cert['radius_upper'])==RU and cert['lambda_scale']==QL and cert['flow_scale']==DF
 assert cert['residual_sha256']==sha(ENUM/'metric_residuals.json')
 check_residual_partition([case['residual_index'] for case in cert['cases']],len(rr))
 covered=set();total=0;counts=Counter();mindual=None;maxdepth=0
 for case in cert['cases']:
  k=case['residual_index'];assert type(k)is int and 0<=k<len(rr) and k!=685 and k not in covered;covered.add(k)
  for key in ('B','I','idx','faces'):assert case[key]==rr[k][key]
  r=replay_tree(case);total+=r['nodes'];counts.update(r['counts']);maxdepth=max(maxdepth,r['maxdepth'])
  if r['mindual'] is not None:mindual=r['mindual'] if mindual is None else min(mindual,r['mindual'])
 assert covered==set(range(len(rr)))-{685}
 assert counts['LOCAL']==0 and mindual is not None and mindual>0
 out={'verified':True,'cases':len(covered),'nodes':total,'counts':dict(counts),'maxdepth':maxdepth,'minimum_radius_margin':str(mindual),'minimum_radius_margin_float':float(mindual),'certificate_sha256':sha(path),'residual_sha256':sha(ENUM/'metric_residuals.json'),'seconds':time.time()-st}
 (BASE/'noncandidate_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True);return out

def verify_candidate():
 st=time.time();rr=json.loads((ENUM/'metric_residuals.json').read_text());path=BASE/'candidate_certificate.json.gz'
 with gzip.open(path,'rt') as f:case=json.load(f)
 assert case['complete'] is True and case['residual_index']==685
 assert rr[685]['B']==case['B']==10 and rr[685]['I']==case['I']==4 and case['faces']==rr[685]['faces']==[list(f) for f in ALL_FACES]
 r=replay_tree(case,True);assert r['counts']['LOCAL']>0 and r['counts']['DUAL']>0
 out={'verified':True,'residual_index':685,'enumeration_idx':rr[685]['idx'],'nodes':r['nodes'],'counts':dict(r['counts']),'maxdepth':r['maxdepth'],'max_local_anchor_distance_squared':str(r['maxlocal']),'local_radius_squared':str(RADIUS**2),'min_dual_radius_margin':str(r['mindual']),'min_dual_radius_margin_float':float(r['mindual']),'root_sha256':sha(BASE/'root_certificate.json'),'certificate_sha256':sha(path),'seconds':time.time()-st}
 (BASE/'candidate_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True);return out
if __name__=='__main__':
 import sys
 if len(sys.argv)>1 and sys.argv[1]=='candidate':verify_candidate()
 else:verify_noncandidate()
