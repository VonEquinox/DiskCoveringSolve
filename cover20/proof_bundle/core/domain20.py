"""Necessary angular domains; rebuilt from every full center-face set."""
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Run without -O/PYTHONOPTIMIZE')
N=20;V=58;D=10000000;A=873043;C4=1822179;C6=3019184
RU=F('0.27084817843999')
def edges(B,faces):
 assert type(B)is int and 12<=B<20 and len(faces)==38-B
 assert len(set(tuple(f) for f in faces))==len(faces)
 assert all(len(f)==3 and list(f)==sorted(set(f)) and all(type(v)is int and 0<=v<N for v in f) for f in faces)
 return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+k,B+c) for k,fa in enumerate(faces) for c in fa]
def bounds(B,faces):
 edges(B,faces)
 adj=[set() for _ in range(N)]
 for f in faces:
  for u in f:adj[u].update(v for v in f if v!=u)
 mat=[[10**9]*B for _ in range(B)]
 for i in range(B):
  mat[i][i]=0;j=(i+1)%B;mat[i][j]=min(mat[i][j],A-(D if j==0 else 0));mat[j][i]=min(mat[j][i],D if j==0 else 0)
 for i in range(B):
  reach=[{i,(i+1)%B}]
  for _ in range(2):reach.append(reach[-1]|set().union(*(adj[v] for v in reach[-1])))
  for j in range(i+1,B):
   ends={j,(j+1)%B};h=next((h for h in range(3) if reach[h]&ends),None)
   if h is None:continue
   cap=[A,C4,C6][h];k=j-i
   if min(k,B-k)*A<=cap:continue
   fw=k*A<D-cap;bw=(B-k)*A<D-cap
   assert not(fw and bw),'contradictory necessary short arcs'
   if fw:mat[i][j]=min(mat[i][j],cap)
   if bw:mat[j][i]=min(mat[j][i],cap-D)
 for k in range(B):
  for i in range(B):
   for j in range(B):mat[i][j]=min(mat[i][j],mat[i][k]+mat[k][j])
 assert all(mat[i][i]==0 for i in range(B)),'negative necessary cycle'
 return mat
def rootbox(mat):
 return [F(-mat[i][0],D) for i in range(1,len(mat))],[F(mat[0][i],D) for i in range(1,len(mat))]
def tighten(mat,lo,hi):
 # Put rational endpoints on an exact common integer lattice. This is an
 # equality-preserving optimization, not fixed-precision rounding.
 from math import lcm
 B=len(mat);assert len(lo)==len(hi)==B-1
 l=[F(0)]+list(map(F,lo));h=[F(0)]+list(map(F,hi));den=D
 for x in l+h:den=lcm(den,x.denominator)
 li=[x.numerator*(den//x.denominator) for x in l]
 ui=[x.numerator*(den//x.denominator) for x in h];scale=den//D
 lower=[max(li[i]-mat[j][i]*scale for i in range(B)) for j in range(B)]
 upper=[min(ui[i]+mat[i][j]*scale for i in range(B)) for j in range(B)]
 if any(x>y for x,y in zip(lower,upper)):return None
 assert lower[0]==upper[0]==0
 return [F(x,den) for x in lower[1:]],[F(x,den) for x in upper[1:]]
