#!/usr/bin/env python3
"""Exact global convexity certificate for the D3 candidate topology."""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time,pickle
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from exact_core import sincos_interval,PI_L
DELTA=F(73891,100000)
PREC=220;QFIX=1<<PREC

def floordiv(a,b):return a//b
def ceildiv(a,b):return -((-a)//b)
class Iv:
 __slots__=('l','u')
 def __init__(self,l,u=None):self.l=int(l);self.u=int(l if u is None else u);assert self.l<=self.u
 @staticmethod
 def from_frac(a,b=None):
  if b is None:b=a
  a=F(a);b=F(b);return Iv((a*QFIX).numerator//(a*QFIX).denominator,ceildiv((b*QFIX).numerator,(b*QFIX).denominator))
 def __add__(self,o):o=iv(o);return Iv(self.l+o.l,self.u+o.u)
 __radd__=__add__
 def __neg__(self):return Iv(-self.u,-self.l)
 def __sub__(self,o):return self+(-iv(o))
 def __rsub__(self,o):return iv(o)-self
 def __mul__(self,o):
  o=iv(o);p=[self.l*o.l,self.l*o.u,self.u*o.l,self.u*o.u];return Iv(min(p)//QFIX,ceildiv(max(p),QFIX))
 __rmul__=__mul__
 def inv(self):
  assert self.l>0
  return Iv((QFIX*QFIX)//self.u,ceildiv(QFIX*QFIX,self.l))
 def __truediv__(self,o):return self*iv(o).inv()
 def lo(self):return F(self.l,QFIX)
 def hi(self):return F(self.u,QFIX)
 def width(self):return F(self.u-self.l,QFIX)
 def __repr__(self):return f'[{self.l/QFIX:.9g},{self.u/QFIX:.9g}]'
def iv(x):return x if isinstance(x,Iv) else Iv.from_frac(F(x))
ZERO=Iv(0)

def load_box():
 d=json.load(open(ROOT/'cover12_symmetric_kkt_cert.json'));Q=int(d['Qx']);rho=F(int(d['rho_num']),int(d['rho_den']));m={n:F(int(q),Q) for n,q in zip(d['names'],d['xnum'])};X={n:Iv.from_frac(q-rho,q+rho) for n,q in m.items()};return X,m,rho

def add(adj,a,b,w):
 if a==b:return
 adj[a][b]=adj[a].get(b,ZERO)+w;adj[b][a]=adj[b].get(a,ZERO)+w

def build_kron(W):
 A=[f'q{i}' for i in range(9)];C=[f'c{i}' for i in range(12)];P=[f'p{i}' for i in range(6)];adj={v:{} for v in A+C+P}
 pat=[('P','Q'),('S','S'),('Q','P')]*3
 for i,(x,y) in enumerate(pat):add(adj,f'q{i}',f'c{i}',W[x]);add(adj,f'q{i}',f'c{(i+1)%9}',W[y])
 faces=[(1,2,11),(4,5,9),(7,8,10),(0,10,11),(3,9,11),(6,9,10)]
 for k,fa in enumerate(faces):
  for c in fa:
   typ=('B' if c>=9 else 'A') if k<3 else ('D' if c>=9 else 'C');add(adj,f'p{k}',f'c{c}',W[typ])
 for v in P+C:
  ne=list(adj[v].items());S=sum((z for _,z in ne),ZERO);assert S.l>0
  for i,(a,x) in enumerate(ne):
   for b,y in ne[i+1:]:add(adj,a,b,x*y/S)
  for a,_ in ne:adj[a].pop(v,None)
  del adj[v]
 return {(i,j):adj[f'q{i}'].get(f'q{j}',ZERO) for i in range(9) for j in range(i+1,9)}

def ldl(A):
 n=len(A);L=[[F(0)]*n for _ in range(n)];D=[F(0)]*n
 for i in range(n):
  L[i][i]=1
  for j in range(i):L[i][j]=(A[i][j]-sum((L[i][k]*D[k]*L[j][k] for k in range(j)),F(0)))/D[j]
  D[i]=A[i][i]-sum((L[i][k]*L[i][k]*D[k] for k in range(i)),F(0));assert D[i]>0,(i,float(D[i]))
 return D

def verify():
 t=time.time();X,m,rho=load_box();assert 8*DELTA<2*PI_L and 4*DELTA<PI_L
 # The hard-coded D3 graph is exactly the sole candidate residual orbit.
 full=[(0,1,11),(0,8,10),(0,10,11),(1,2,11),(2,3,11),(3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),(6,9,10),(7,8,10),(9,10,11)]
 from enumerate_orbits import canonical_orbit,canon_state
 reps=pickle.load(open(ROOT/'triangulations_B9_I3_orbits.pkl','rb'));assert canonical_orbit(canon_state(full),9,3)==reps[21015]
 metric=json.load(open(ROOT/'metric_exact_cert_safe.json'));assert 21015 in metric['cases']['9,3']['residual_indices']
 # Exact root branch and candidate gaps: h and k=pi/3-2h lie in (0,pi/2),
 # and 2h,2k are both below DELTA.
 c,s,z=X['c'],X['s'],X['z'];half=iv(F(1,2));L=(c*c-s*s)*half+z*c*s;M=z*(c*c-s*s)*half-c*s
 sinDl=sincos_interval(DELTA/2)[0];assert min(c.l,s.l,L.l,M.l)>0 and s.hi()<sinDl and M.hi()<sinDl
 W={k:X[k] for k in ['P','Q','S','A','B','C','D']};assert min(z.l for z in W.values())>0
 K=build_kron(W);assert min(z.l for z in K.values())>=0
 coslo={d:sincos_interval(d*DELTA)[2] for d in range(1,5)}
 assert coslo[1]>0 and coslo[2]>0 and coslo[3]<0 and coslo[4]<0
 H=[[F(0)]*9 for _ in range(9)]
 for (i,j),kap in K.items():
  d=min(j-i,9-j+i);cl=coslo[d];a=2*(kap.lo() if cl>=0 else kap.hi())*cl
  H[i][i]+=a;H[j][j]+=a;H[i][j]-=a;H[j][i]-=a
 piv=ldl([row[1:] for row in H[1:]])
 report={'verified':True,'fixed_precision_bits':PREC,'delta':str(DELTA),'cos_lower':{str(k):str(v) for k,v in coslo.items()},'min_conductance_lower':float(min(z.lo() for z in K.values())),'max_conductance_width':float(max(z.width() for z in K.values())),'min_ldl_pivot':float(min(piv)),'ldl_pivots':[float(x) for x in piv],'candidate_orbit':21015,'gap_half_sine_margin':float(min(sinDl-s.hi(),sinDl-M.hi())),'seconds':time.time()-t}
 json.dump(report,open(ROOT/'cover12_candidate_convex_verified.json','w'),indent=2);print(json.dumps(report,indent=2));return report
if __name__=='__main__':verify()
