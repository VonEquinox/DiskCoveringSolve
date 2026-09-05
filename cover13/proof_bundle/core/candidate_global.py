#!/usr/bin/env python3
"""Generate/replay a finite, exact optimality certificate for ONE specified
13-center triangulation. This is NOT an all-topologies covering theorem.
Only Python standard-library exact arithmetic is used.
"""
from pathlib import Path
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
from fractions import Fraction as F
from math import isqrt
import json,gzip,hashlib,time,sys
from exact_arcs import SCALE,QA,QB,sincos_turn,arc_halfplane,rational_interval,dot_interval
ROOT=Path(__file__).resolve().parent
CAP=F(1127,10000);QL=10**15;DF=10**36
RU=F(346645456927389644,10**18)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilroot(n):
 assert n>=0;r=isqrt(n);return r+(r*r<n)

def topology():
 root=json.loads((ROOT/'cover13_kkt_root_119d.json').read_text())
 ff=root['all_faces'];assert len(ff)==14 and len(set(tuple(f) for f in ff))==14
 ed=[(i,10+i) for i in range(10)]+[(i,10+(i+1)%10) for i in range(10)]+[(23+k,10+c) for k,f in enumerate(ff) for c in f]
 assert len(ed)==62
 active=root['active_faces'];assert len(active)==9
 original=[tuple(e) for e in root['edges']]
 # Active graph is an exact subgraph of this 14-face graph after witness relabeling.
 amap={23+k:23+ff.index(f) for k,f in enumerate(active)}
 mapped={tuple(sorted((amap.get(u,u),amap.get(v,v)))) for u,v in original}
 assert mapped<=set(tuple(sorted(e)) for e in ed)
 return ff,ed

FACES,EDGES=topology()

def tighten(lo,hi):
 l=[F(0)]+list(lo)+[F(1)];h=[F(0)]+list(hi)+[F(1)]
 for _ in range(12):
  old=(tuple(l),tuple(h))
  for i in range(1,11):h[i]=min(h[i],h[i-1]+CAP);l[i]=max(l[i],l[i-1])
  for i in range(9,-1,-1):h[i]=min(h[i],h[i+1]);l[i]=max(l[i],l[i+1]-CAP)
  if any(a>b for a,b in zip(l,h)):return None
  if old==(tuple(l),tuple(h)):return l[1:10],h[1:10]
 raise AssertionError('unstable interval propagation')

def root_coordinates():
 c=json.loads((ROOT/'cover13_krawczyk_cert.json').read_text());Q=int(c['Qx']);rho=F(int(c['rho_num']),int(c['rho_den']))
 x=[F(int(v),Q) for v in c['xnum']]
 assert RU*RU>x[62]+rho
 # delta=2 asin(r*) < 2*pi*CAP: certify r* < sin(pi*CAP).
 _,sn=sincos_turn(CAP/2)
 assert sn[0]>0 and F(sn[0],SCALE)**2>x[62]+rho
 q=[]
 for i in range(9):
  row=[]
  for j in (2*i,2*i+1):
   lo=rational_interval(x[j]-rho)[0];hi=rational_interval(x[j]+rho)[1];row.append((lo,hi))
  q.append(row)
 return q


def local_check(lo,hi,qs):
 total=0
 for i in range(9):
  assert hi[i]-lo[i]<F(1,2)
  cm,sm=sincos_turn((lo[i]+hi[i])/2)
  # Positive midpoint dot product and angular halfwidth <pi/2 rule out
  # the unique interior minimum (-1). Thus the minimum is at an endpoint.
  assert dot_interval(qs[i],(cm,sm))[0]>0
  vals=[]
  for u in (lo[i],hi[i]):
   c,s=sincos_turn(u);vals.append(dot_interval(qs[i],(c,s))[0])
  total+=2*SCALE-2*min(vals)
 assert 64*total<=SCALE,('local radius failed',float(F(total,SCALE)))
 return F(1,64)-F(total,SCALE)

def dual_check(lo,hi,record):
 lam=record['lambda_num'];flow=record['force_num']
 assert len(lam)==9 and all(type(v) is int and v>=0 for v in lam)
 assert len(flow)==len(EDGES) and all(len(row)==2 and all(type(v) is int for v in row) for row in flow)
 normals=[];bs=[]
 for l,h in zip(lo,hi):
  aa,b=arc_halfplane(l,h);normals.append(aa);bs.append(b)
 target=[[0,0] for _ in range(37)]
 for i in range(9):
  for d in range(2):target[i+1][d]=lam[i]*normals[i][d]*(DF//(QL*QA));target[0][d]-=target[i+1][d]
 actual=[[0,0] for _ in range(37)]
 for (u,v),f in zip(EDGES,flow):
  for d in range(2):actual[u][d]+=f[d];actual[v][d]-=f[d]
 assert actual==target,'force imbalance'
 normsum=sum(ceilroot(f[0]*f[0]+f[1]*f[1]) for f in flow)
 Bnum=sum(lam[i]*(bs[i]-normals[i][0]*(QB//QA)) for i in range(9))
 assert Bnum>0 and normsum>0
 left=Bnum*RU.denominator*DF;right=RU.numerator*normsum*QL*QB
 assert left>right,('nonpositive dual margin',float(F(Bnum,QL*QB)/F(normsum,DF)-RU))
 return F(Bnum*DF,QL*QB*normsum)-RU


def verify():
 st=time.time();qs=root_coordinates();path=ROOT/'candidate_global_certificate.json.gz'
 with gzip.open(path,'rt',encoding='utf8') as f:c=json.load(f)
 assert c['schema']=='candidate-arc-force-tree-v1'
 assert c['scope']=='single-10-boundary-3-interior-triangulation-only'
 assert F(c['cap'])==CAP and c['flow_scale']==DF and c['lambda_scale']==QL and F(c['radius_upper'])==RU
 assert c['all_faces']==FACES
 assert c['root_sha256']==sha(ROOT/'cover13_krawczyk_cert.json')
 assert c['root_data_sha256']==sha(ROOT/'cover13_kkt_root_119d.json')
 assert c['local_certificate_sha256']==sha(ROOT/'anchor_isolation_certificate.json')
 ns=c['nodes'];seen=set();counts={};mdu=None;mloc=None
 lo=[max(F(0),1-(10-i)*CAP) for i in range(1,10)];hi=[min(F(1),i*CAP) for i in range(1,10)]
 stack=[(0,lo,hi,0)];maxdepth=0
 while stack:
  idx,lo,hi,depth=stack.pop();assert type(idx) is int and 0<=idx<len(ns) and idx not in seen
  seen.add(idx);maxdepth=max(maxdepth,depth);node=ns[idx];kind=node['kind'];box=tighten(lo,hi)
  if kind=='EMPTY':assert box is None;counts[kind]=counts.get(kind,0)+1;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   j=node['axis'];mid=F(node['mid']);children=node['children'];assert type(j) is int and 0<=j<9
   assert lo[j]<mid<hi[j] and len(children)==2 and children[0]!=children[1]
   lh=list(hi);lh[j]=mid;rl=list(lo);rl[j]=mid
   stack.extend([(children[1],rl,list(hi),depth+1),(children[0],list(lo),lh,depth+1)])
  elif kind=='LOCAL':
   mar=local_check(lo,hi,qs);mloc=mar if mloc is None else min(mloc,mar);counts[kind]=counts.get(kind,0)+1
  elif kind=='DUAL':
   mar=dual_check(lo,hi,node);mdu=mar if mdu is None else min(mdu,mar);counts[kind]=counts.get(kind,0)+1
  else:raise AssertionError(('unknown node kind',kind))
 assert seen==set(range(len(ns)))
 out={'verified':True,'scope':c['scope'],'nodes':len(ns),'leaves':counts,'maxdepth':maxdepth,
 'minimum_dual_radius_margin':str(mdu),'minimum_dual_radius_margin_float':float(mdu),
 'minimum_local_squared_radius_margin':str(mloc),'minimum_local_squared_radius_margin_float':float(mloc),
 'radius_upper':str(RU),'certificate_sha256':sha(path),'seconds':time.time()-st,
 'standalone_scope':'candidate topology only; the master verifier covers every other topology'}
 (ROOT/'candidate_global_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 verify()
