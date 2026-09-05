"""Exact difference-bound domains implied by boundary and four-rod arc bounds."""
from fractions import Fraction as F
from verify_partition import A,C,D,signature

def bounds(B,faces):
 M=[[10**9]*B for _ in range(B)]
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
