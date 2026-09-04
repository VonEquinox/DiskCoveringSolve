#!/usr/bin/env python3
from fractions import Fraction as F
from collections import deque,defaultdict
from pathlib import Path
import itertools,pickle,json,math,time,sys,hashlib
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parent))
from exact_core import PI_L,PI_U,sincos_interval
ROOT=Path(__file__).resolve().parent
T=F(13039536,10**8)
ALPHA={2:F(73891,100000),3:F(57241,50000),4:F(161399,100000),5:F(22523,10000)}
Q=10**12
CASES=[(9,1),(9,2),(9,3),(10,1),(10,2),(11,1)]

def canon_face(f):return tuple(sorted(map(int,f)))
def canon_state(fs):return tuple(sorted(canon_face(f) for f in fs))
def transform_state(s,B,I,shift,reflect,perm):
 def mp(x):return ((-x if reflect else x)+shift)%B if x<B else B+perm[x-B]
 return canon_state(tuple(mp(x) for x in f) for f in s)
def orbit_set(s,B,I):
 return {transform_state(s,B,I,sh,refl,p) for refl in (False,True) for sh in range(B) for p in (itertools.permutations(range(I)) if I else [()])}
def valid_disk_triangulation(s,B,I):
 s=canon_state(s);V=B+I;Fcnt=B+2*I-2;Ecnt=2*B+3*I-3
 if len(s)!=Fcnt or len(set(s))!=len(s):return False
 if any(len(set(f))<3 or min(f)<0 or max(f)>=V for f in s):return False
 if set(x for f in s for x in f)!=set(range(V)):return False
 boundary={tuple(sorted((i,(i+1)%B))) for i in range(B)}
 ec=defaultdict(int);ef=defaultdict(list)
 for fi,(a,b,c) in enumerate(s):
  for x,y in ((a,b),(a,c),(b,c)):
   e=tuple(sorted((x,y)));ec[e]+=1;ef[e].append(fi)
 if len(ec)!=Ecnt:return False
 for e,n in ec.items():
  if e in boundary:
   if n!=1:return False
  elif n!=2:return False
 if any(ec.get(e,0)!=1 for e in boundary):return False
 # Face adjacency through interior edges must be connected.
 adj=[set() for _ in s]
 for e,ff in ef.items():
  if len(ff)==2:adj[ff[0]].add(ff[1]);adj[ff[1]].add(ff[0])
 seen={0};stack=[0]
 while stack:
  q=stack.pop()
  for z in adj[q]:
   if z not in seen:seen.add(z);stack.append(z)
 if len(seen)!=len(s):return False
 # Every vertex link is a single circle (interior vertex) or a single path
 # with the prescribed boundary neighbours as endpoints (boundary vertex).
 for x in range(V):
  ladj=defaultdict(set);ledges=set()
  for f in s:
   if x not in f:continue
   a,b=[q for q in f if q!=x];e=tuple(sorted((a,b)))
   if e in ledges:return False
   ledges.add(e);ladj[a].add(b);ladj[b].add(a)
  if not ladj:return False
  root=next(iter(ladj));vis={root};st=[root]
  while st:
   q=st.pop()
   for z in ladj[q]:
    if z not in vis:vis.add(z);st.append(z)
  if len(vis)!=len(ladj):return False
  deg={q:len(z) for q,z in ladj.items()}
  if x<B:
   ends={q for q,d in deg.items() if d==1}
   if ends!={(x-1)%B,(x+1)%B}:return False
   if any(d not in (1,2) for d in deg.values()):return False
  elif any(d!=2 for d in deg.values()):return False
 # A connected compact 2-manifold with this one boundary cycle and
 # V-E+F=1 is a topological disk.
 return V-len(ec)+len(s)==1

def brown_labeled(B,I):
 m=B-3;n=I
 z=F(2*math.factorial(2*m+3)*math.factorial(4*n+2*m+1),math.factorial(m+2)*math.factorial(m)*math.factorial(n)*math.factorial(3*n+2*m+3))
 z*=math.factorial(I);assert z.denominator==1;return z.numerator

def graph(B,I,faces):
 ncent=B+I;n=B+ncent+len(faces);adj=[[] for _ in range(n)]
 def add(a,b):adj[a].append(b);adj[b].append(a)
 for i in range(B):add(i,B+i);add(i,B+(i+1)%B)
 for k,f in enumerate(faces):
  p=B+ncent+k
  for v in f:add(p,B+v)
 return adj

def anchor_dist(B,I,faces):
 adj=graph(B,I,faces);out=[]
 for s in range(B):
  dd=[99]*len(adj);dd[s]=0;q=deque([s])
  while q:
   u=q.popleft()
   for v in adj[u]:
    if dd[v]>dd[u]+1:dd[v]=dd[u]+1;q.append(v)
  out.append(tuple(dd[:B]))
 return tuple(out)
def keyD(D):return ','.join(''.join(str(x) for x in row) for row in D)
def constraints(D):
 B=len(D);rows=[];rhs=[]
 for i in range(B):
  for j in range(i+1,B):
   k=j-i
   if 2*k==B:continue
   inds=list(range(i,j)) if 2*k<B else list(range(j,B))+list(range(i))
   # The chosen combinatorial path is certified to be the minor circular arc only when its cap is < pi.
   if len(inds)*ALPHA[2]>=PI_L:continue
   d=D[i][j]
   if d<=5:
    row=[0]*B
    for h in inds:row[h]=1
    rows.append(row);rhs.append(ALPHA[d])
 return rows,rhs

def make_cert(D):
 rows,rhs=constraints(D);B=len(D);A=np.array(rows,float);b=np.array([float(x) for x in rhs])
 if not len(rows):return None
 res=linprog(b,A_ub=-A.T,b_ub=-np.ones(B),bounds=[(0,None)]*len(rows),method='highs')
 if not res.success or res.fun>=2*math.pi-1e-8:return None
 nums=[max(0,math.ceil(float(x)*Q+1e-7)) for x in res.x]
 cov=[sum(nums[j]*rows[j][i] for j in range(len(rows))) for i in range(B)]
 for i in range(B):
  if cov[i]<Q:
   j=next(j for j,row in enumerate(rows) if row[i]);d=Q-cov[i];nums[j]+=d
   for h in range(B):cov[h]+=d*rows[j][h]
 cost=sum((F(n,Q)*rhs[j] for j,n in enumerate(nums)),F(0))
 if cost>=2*PI_L:return None
 return {'lam':[[j,n] for j,n in enumerate(nums) if n],'cost':[cost.numerator,cost.denominator]}

def verify_angle():
 out={}
 for d,a in ALPHA.items():
  sl,su,cl,cu=sincos_interval(a/2);mar=4*sl*sl-d*d*T;assert mar>0;out[str(d)]=[mar.numerator,mar.denominator]
 assert 8*ALPHA[2]<2*PI_L
 assert 4*ALPHA[2]<PI_L
 return out

def generate(path):
 angle=verify_angle();data={'version':1,'T':[T.numerator,T.denominator],'alpha':{str(k):[v.numerator,v.denominator] for k,v in ALPHA.items()},'angle_margin':angle,'cases':{}}
 for B,I in CASES:
  reps=pickle.load(open(ROOT/f'triangulations_B{B}_I{I}_orbits.pkl','rb'));cache={};res=[];t=time.time()
  for idx,s in enumerate(reps):
   D=anchor_dist(B,I,s);k=keyD(D)
   if k not in cache:cache[k]=make_cert(D)
   if cache[k] is None:res.append(idx)
  data['cases'][f'{B},{I}']={'orbit_count':len(reps),'certs':{k:v for k,v in cache.items() if v is not None},'residual_indices':res}
  print('GEN',(B,I),'orbits',len(reps),'uniqueD',len(cache),'residual',len(res),'sec',time.time()-t,flush=True)
 json.dump(data,open(path,'w'),separators=(',',':'));print('saved',path)

def verify(path):
 t0=time.time();angle=verify_angle();data=json.load(open(path));assert F(*data['T'])==T;total=0;resall=0;reports=[]
 for B,I in CASES:
  reps=pickle.load(open(ROOT/f'triangulations_B{B}_I{I}_orbits.pkl','rb'));assert reps==sorted(set(reps));seen_count=0
  for s in reps:
   assert valid_disk_triangulation(s,B,I);orb=orbit_set(s,B,I);assert min(orb)==s;seen_count+=len(orb)
  assert seen_count==brown_labeled(B,I),(B,I,seen_count,brown_labeled(B,I))
  blk=data['cases'][f'{B},{I}'];assert blk['orbit_count']==len(reps);res=set(map(int,blk['residual_indices']));ne=0;minmar=None
  for idx,s in enumerate(reps):
   D=anchor_dist(B,I,s);k=keyD(D);rows,rhs=constraints(D)
   if k in blk['certs']:
    z=blk['certs'][k];nums=[0]*len(rows)
    for j,n in z['lam']:nums[int(j)]=int(n)
    cov=[sum(nums[j]*rows[j][i] for j in range(len(rows))) for i in range(B)];assert min(cov)>=Q
    cost=sum((F(nums[j],Q)*rhs[j] for j in range(len(rows))),F(0));assert cost==F(*z['cost']);mar=2*PI_L-cost;assert mar>0
    minmar=mar if minmar is None or mar<minmar else minmar;ne+=1
   else:assert idx in res
  assert ne+len(res)==len(reps);total+=len(reps);resall+=len(res)
  reports.append({'B':B,'I':I,'orbits':len(reps),'labeled':seen_count,'excluded':ne,'residual':len(res),'min_margin':float(minmar)})
  print('VER',(B,I),reports[-1],flush=True)
 out={'verified':True,'T':str(T),'certificate':str(Path(path).name),'total_orbits':total,'total_residual':resall,'cases':reports,'seconds':time.time()-t0};report_name='metric_exact_safe_verified.json' if 'safe' in Path(path).stem else 'metric_exact_verified.json';json.dump(out,open(ROOT/report_name,'w'),indent=2);print(json.dumps(out,indent=2))
if __name__=='__main__':
 mode=sys.argv[1];path=sys.argv[2] if len(sys.argv)>2 else str(ROOT/'metric_exact_cert.json')
 generate(path) if mode=='generate' else verify(path)
