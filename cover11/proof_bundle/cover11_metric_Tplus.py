#!/usr/bin/env python3
"""Exact corrected shortest-path/Farkas screen at T+=0.1443877286317993."""
from fractions import Fraction as F
from collections import deque
import pickle,json,math,sys,numpy as np
from scipy.optimize import linprog
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core
from cover11_fixed_trig import Q as TQ,enc,sin_iv
T=core.TARGET
ALPHA={2:F(1949,2500),3:F(1213,1000),4:F(8633,5000),5:F(125311,50000)}
Q=10**10
CASES=[(9,1,'/mnt/data/triangulations_B9_I1.pkl'),(9,2,'/mnt/data/triangulations_B9_I2.pkl'),(10,1,'/mnt/data/triangulations_B10_I1.pkl')]
def graph(B,I,faces):
 ncent=B+I;n=B+ncent+len(faces);edges=[]
 for i in range(B):edges += [(i,B+i),(i,B+(i+1)%B)]
 for k,f in enumerate(faces):
  p=B+ncent+k;edges += [(p,B+v) for v in f]
 return n,edges
def anchor_dist(B,I,faces):
 n,edges=graph(B,I,faces);adj=[[] for _ in range(n)]
 for a,b in edges:adj[a].append(b);adj[b].append(a)
 D=[[0]*B for _ in range(B)]
 for s in range(B):
  dd=[999]*n;dd[s]=0;q=deque([s])
  while q:
   u=q.popleft()
   for v in adj[u]:
    if dd[v]>dd[u]+1:dd[v]=dd[u]+1;q.append(v)
  D[s]=dd[:B]
 return D
def keyD(D):return ''.join(chr(48+x) for row in D for x in row)
def constraints(D):
 B=len(D);rows=[];rhs=[]
 for i in range(B):
  for j in range(i+1,B):
   k=j-i
   if 2*k==B:continue
   inds=list(range(i,j)) if 2*k<B else list(range(j,B))+list(range(i))
   d=D[i][j]
   if d<=5:
    row=[0]*B
    for h in inds:row[h]=1
    rows.append(row);rhs.append(ALPHA[d])
 return rows,rhs
def make_cert(D):
 rows,rhs=constraints(D);B=len(D);A=np.array(rows,float);b=np.array([float(x) for x in rhs])
 res=linprog(b,A_ub=-A.T,b_ub=-np.ones(B),bounds=[(0,None)]*len(rows),method='highs-ds',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
 assert res.success
 if res.fun>=2*math.pi-1e-11:return None
 nums=[max(0,math.ceil(float(x)*Q+1e-6)) for x in res.x]
 cov=[sum(nums[j]*rows[j][i] for j in range(len(rows))) for i in range(B)]
 for i in range(B):
  if cov[i]<Q:
   j=next(j for j,row in enumerate(rows) if row[i]);add=Q-cov[i];nums[j]+=add
   for h in range(B):cov[h]+=add*rows[j][h]
 assert min(cov)>=Q
 cost=sum(F(n,Q)*rhs[j] for j,n in enumerate(nums))
 if not cost<core.S_L:
  raise RuntimeError(('rounded Farkas cost',float(cost),float(core.S_L),float(res.fun)))
 return {'q':Q,'lam':[[j,n] for j,n in enumerate(nums) if n],'cost':[cost.numerator,cost.denominator]}
def angle_verify():
 for d,a in ALPHA.items():
  iv=sin_iv(enc(a/2));sl=F(iv.lo,TQ)
  assert sl>0 and 4*sl*sl>d*d*T,(d,float(4*sl*sl-d*d*T))
 assert 4*ALPHA[2]<core.PI_L

def generate(path):
 angle_verify();data={'version':2,'target':[T.numerator,T.denominator],'alpha':{str(d):[x.numerator,x.denominator] for d,x in ALPHA.items()},'cases':{}}
 for B,I,fn in CASES:
  reps=pickle.load(open(fn,'rb'));certs={};resid=[]
  for idx,faces in enumerate(reps):
   D=anchor_dist(B,I,faces);key=keyD(D)
   if key not in certs:certs[key]=make_cert(D)
   if certs[key] is None:resid.append({'idx':idx,'faces':[list(f) for f in faces],'D':key})
  exc={k:v for k,v in certs.items() if v is not None}
  data['cases'][f'{B},{I}']={'orbit_count':len(reps),'certs':exc,'residual':resid}
  print((B,I),'orbits',len(reps),'uniqueD',len(certs),'excludedD',len(exc),'residual',len(resid),flush=True)
 json.dump(data,open(path,'w'),separators=(',',':'))

def verify(path):
 angle_verify();data=json.load(open(path));assert F(*data['target'])==T
 for d,a in ALPHA.items():assert F(*data['alpha'][str(d)])==a
 for B,I,fn in CASES:
  reps=pickle.load(open(fn,'rb'));block=data['cases'][f'{B},{I}'];assert len(reps)==int(block['orbit_count'])
  residual={(int(z['idx']),tuple(tuple(f) for f in z['faces'])) for z in block['residual']};ne=nr=0;mn=None
  for idx,faces in enumerate(reps):
   D=anchor_dist(B,I,faces);key=keyD(D);rows,rhs=constraints(D)
   if key in block['certs']:
    z=block['certs'][key];q=int(z['q']);nums=[0]*len(rows)
    for j,n in z['lam']:nums[int(j)]=int(n)
    cov=[sum(nums[j]*rows[j][h] for j in range(len(rows))) for h in range(B)];assert min(cov)>=q
    cost=sum(F(nums[j],q)*rhs[j] for j in range(len(rows)));assert cost==F(*z['cost']) and cost<core.S_L
    mar=core.S_L-cost;mn=mar if mn is None or mar<mn else mn;ne+=1
   else:
    assert (idx,tuple(faces)) in residual;nr+=1
  assert nr==len(residual)
  print('VERIFIED',(B,I),'excluded',ne,'residual',nr,'minmargin',float(mn),flush=True)
if __name__=='__main__':
 mode=sys.argv[1];path=sys.argv[2] if len(sys.argv)>2 else '/mnt/data/cover11_metric_Tplus.json'
 generate(path) if mode=='generate' else verify(path)
