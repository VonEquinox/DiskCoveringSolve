"""Certify exact D2 symmetry using equivariance and uniqueness of the KKT root.
The interval image of the tight root box is checked inside the large uniqueness box.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
import system14 as S
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[a0*b0 for a0 in a for b0 in b];return min(v),max(v)
def point(x):x=F(x);return x,x
def scale(a,x):return mul(a,point(x))
def vadd(a,b):return [add(x,y) for x,y in zip(a,b)]
def vsub(a,b):return [sub(x,y) for x,y in zip(a,b)]
def dot(a,b):return add(mul(a[0],b[0]),mul(a[1],b[1]))
def cross(a,b):return sub(mul(a[0],b[1]),mul(a[1],b[0]))
def norm2(a):return dot(a,a)

def rootbox(rho=F(1,10**80)):
 d=json.load(open(R/'root14.json'));q=int(d['Qx']);z=[F(int(v),q) for v in d['xnum']];iv=[(v-rho,v+rho) for v in z];p=[[point(1),point(0)]]+[iv[2*i:2*i+2] for i in range(S.G//2)]
 return d,z,iv,p

def automorphism(kind):
 if kind=='halfturn':qc=[(i+5)%10 for i in range(10)];cc=qc+[12,13,10,11];fixed_source=5
 elif kind=='reflection':qc=[(-i-1)%10 for i in range(10)];cc=[(-i)%10 for i in range(10)]+[10,13,12,11];fixed_source=9
 else:raise ValueError(kind)
 ad={tuple(sorted(f)):i for i,f in enumerate(S.ACTIVE)}
 pp=[ad[tuple(sorted(cc[v] for v in f))] for f in S.ACTIVE]
 perm=qc+[10+c for c in cc]+[24+p for p in pp]
 assert sorted(perm)==list(range(32)) and all(perm[perm[i]]==i for i in range(32)) and perm[fixed_source]==0
 ed={tuple(sorted(e)):i for i,e in enumerate(S.EDGES)};ep=[ed[tuple(sorted((perm[u],perm[v])))] for u,v in S.EDGES]
 assert sorted(ep)==list(range(S.M))
 return perm,ep,fixed_source

def verify():
 st=time.time();d,z,iv,p=rootbox();mu=[None]+iv[S.P+S.M:];weights=iv[S.P:S.P+S.M]
 # Full unpinned stationarity at q0 follows from rotational invariance. This is
 # its uniquely determined radial multiplier: mu0=-sum w(1-neighbor_x).
 m0=point(0)
 for w,(u,v) in zip(weights,S.EDGES):
  if u==0:m0=sub(m0,mul(w,sub(point(1),p[v][0])))
  elif v==0:m0=sub(m0,mul(w,sub(point(1),p[u][0])))
 mu[0]=m0
 reports=[]
 for kind in ['halfturn','reflection']:
  perm,ep,j=automorphism(kind);a,b=p[j]
  T=[[a,b],[neg(b),a]] if kind=='halfturn' else [[a,b],[b,neg(a)]]
  image=[None]*S.D;geo=[None]*32
  for old,new in enumerate(perm):
   v=[dot(row,p[old]) for row in T];geo[new]=v
   if new:
    for k in range(2):image[2*(new-1)+k]=v[k]
  image[S.TID]=iv[S.TID]
  for i,jj in enumerate(ep):image[S.P+jj]=weights[i]
  for i in range(10):
   j=perm[i]
   if j:image[S.P+S.M+j-1]=mu[i]
  assert all(v is not None for v in image)
  deviations=[max(abs(v[0]-c),abs(v[1]-c)) for v,c in zip(image,z)]
  assert max(deviations)<F(1,10**60)
  reports.append({'generator':kind,'node_permutation':perm,'edge_permutation':ep,'image_max_midpoint_distance':str(max(deviations))})
 # The half-turn involution is not the identity, hence it is exactly -I.
 assert p[5][0][1]<0 and p[0][0][0]>0
 # Reflection fixes the c0/c10/c12 axis; c0 is nonzero.
 assert norm2(p[10])[0]>F(1,2)
 # The two degree-two inner vertices have equal positive rod weights under
 # reflection, as does c0. Their transverse offsets are equal to +/- n.
 perm,ep,j=automorphism('reflection')
 for center,wits in [(10,[0,9]),(20,[25,31]),(22,[27,29])]:
  inc=[(k,v if u==center else u) for k,(u,v) in enumerate(S.EDGES) if u==center or v==center]
  assert sorted(v for k,v in inc)==sorted(wits) and len(inc)==2
  assert ep[inc[0][0]]==inc[1][0] and perm[center]==center
 n=vsub(p[0],p[10]);n1=vsub(p[25],p[20]);n3=vsub(p[27],p[22])
 assert dot(n,n1)[0]>0 and dot(n,n3)[0]>0
 out={'verified':True,'generators':reports,'root_sha256':hashlib.sha256((R/'root14.json').read_bytes()).hexdigest(),
      'consequences':['halfturn is -I','c0=(q0+q9)/2','p1-c10=p3-c12=q0-c0','p7-c10=p5-c12=-(q0-c0)'],'seconds':time.time()-st}
 (R/'symmetry14_verified.json').write_text(json.dumps(out,indent=2));print('SYMMETRY VERIFIED',float(max(F(r['image_max_midpoint_distance']) for r in reports)),flush=True)
 return out
if __name__=='__main__':verify()
