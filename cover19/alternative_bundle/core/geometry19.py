"""Explicit Q(sqrt(3)) candidate; stored pair (x,y) represents (x,sqrt(3)*y).
No optimizer or algebraic-root assertion is trusted. All geometry and weights
are rational after the invertible Cartesian scaling diag(1,sqrt(3)).
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from math import isqrt
if not __debug__:raise RuntimeError('Exact checks require assertions enabled')
R=Path(__file__).resolve().parent
B=12;N=19;T=F(1,13)
LATTICE=[(2,0),(1,1),(0,2),(-1,2),(-2,2),(-2,1),(-2,0),(-1,-1),(0,-2),(1,-2),(2,-2),(2,-1),(0,0),(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def scale(a,t):return (t*a[0],t*a[1])
def dot(a,b):return a[0]*b[0]+3*a[1]*b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0] # physical determinant / sqrt(3)
def norm2(a):return dot(a,a)
CENTERS=[(F(12*m+9*n,26),F(-2*m+5*n,26)) for m,n in LATTICE]
# Q_0=(1,0). For adjacent centers a,b the outside intersection of the
# radius-1/sqrt(13) circles is their midpoint plus J_clockwise(b-a)/(2 sqrt(3)).
# In the stored pair coordinates this offset is (d_y/2,-d_x/6).
ANCHORS=[]
for i in range(B):
 a,b=CENTERS[i],CENTERS[(i+1)%B];d=sub(b,a)
 ANCHORS.append(add(scale(add(a,b),F(1,2)),(d[1]/2,-d[0]/6)))
# All elementary equilateral triangular lattice faces, computed from centers.
FACES=[f for f in combinations(range(N),3) if all(norm2(sub(CENTERS[i],CENTERS[j]))==3*T for i,j in combinations(f,2))]
ACTIVE=[(1,13,14),(3,14,15),(5,15,16),(7,16,17),(9,17,18),(11,13,18)]
FACEPOINTS=[scale(add(add(CENTERS[a],CENTERS[b]),CENTERS[c]),F(1,3)) for a,b,c in FACES]
# The center at the origin is deliberately omitted from the stressed subgraph.
# Its location is irrelevant to the lower-bound implication, not to the upper cover.
CENTER_IDS=[i for i in range(N) if i!=12]
CENTER_MAP={c:B+i for i,c in enumerate(CENTER_IDS)}
POSITIVE_POINTS=[scale(add(add(CENTERS[a],CENTERS[b]),CENTERS[c]),F(1,3)) for a,b,c in ACTIVE]
Z=ANCHORS+[CENTERS[i] for i in CENTER_IDS]+POSITIVE_POINTS
EDGES=[(i,CENTER_MAP[i]) for i in range(B)]+[(i,CENTER_MAP[(i+1)%B]) for i in range(B)]+[(B+len(CENTER_IDS)+k,CENTER_MAP[c]) for k,fa in enumerate(ACTIVE) for c in fa]
WEIGHTS=[F(1,52) if i%2==0 else F(1,39) for i in range(B)]+[F(1,39) if i%2==0 else F(1,52) for i in range(B)]+[F(1,39)]*18
MUS=[F(-1,156)]*(B-1)
G=2*(len(Z)-1);M=len(EDGES);C=M+B-1

def check_geometry():
 assert len(set(LATTICE))==N and set(LATTICE)=={(m,n) for m in range(-2,3) for n in range(-2,3) if max(abs(m),abs(n),abs(m+n))<=2}
 assert len(FACES)==24 and len(ACTIVE)==6 and set(ACTIVE)<=set(FACES)
 assert ANCHORS[0]==(F(1),F(0)) and all(norm2(q)==1 for q in ANCHORS)
 assert all(norm2(sub(ANCHORS[i],CENTERS[c]))==T for i in range(B) for c in (i,(i+1)%B))
 assert all(norm2(sub(p,CENTERS[c]))==T for f,p in zip(FACES,FACEPOINTS) for c in f)
 assert G==70 and M==42 and C==53 and len(Z)==36
 assert all(norm2(sub(Z[u],Z[v]))==T for u,v in EDGES)
 assert len(WEIGHTS)==M and len(MUS)==B-1 and min(WEIGHTS)>0 and max(MUS)<0 and sum(WEIGHTS)==1
 A,energy,P=matrices()
 assert all(sum(WEIGHTS[k]*A[k][j] for k in range(M))+sum(MUS[i]*A[M+i][j] for i in range(B-1))==0 for j in range(G))
 return True

def matrices():
 """A is the physical constraint Jacobian times S; energy is S^T M S.
 Here S=diag(1,sqrt(3),...). P is S^T P_anchor S. All entries rational.
 """
 A=[[F(0)]*G for _ in range(C)];energy=[[F(0)]*G for _ in range(G)];P=[F(0)]*G
 for k,(u,v) in enumerate(EDGES):
  for d,fac in ((0,1),(1,3)):
   dif=Z[u][d]-Z[v][d]
   for n,sg in ((u,1),(v,-1)):
    if n:A[k][2*(n-1)+d]+=2*fac*sg*dif
   for n,sg in ((u,1),(v,-1)):
    if n:
     for nn,ss in ((u,1),(v,-1)):
      if nn:energy[2*(n-1)+d][2*(nn-1)+d]+=fac*WEIGHTS[k]*sg*ss
 for i in range(1,B):
  for d,fac in ((0,1),(1,3)):
   A[M+i-1][2*(i-1)+d]=2*fac*Z[i][d]
   energy[2*(i-1)+d][2*(i-1)+d]+=fac*MUS[i-1]
   P[2*(i-1)+d]=F(fac)
 return A,energy,P

def sqrt_interval(n,den=10**70):
 v=isqrt(n*den*den)
 return F(v,den),F(v+1,den)

def cartesian_intervals(points):
 sl,su=sqrt_interval(3)
 out=[]
 for x,y in points:
  yy=(y*sl,y*su) if y>=0 else (y*su,y*sl)
  out.append([(x,x),yy])
 return out
