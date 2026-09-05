"""Independent exact replay primitives. Standard library only."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
from exact_arcs import QA,QB,arc_halfplane
if not __debug__:raise RuntimeError("Run without -O/PYTHONOPTIMIZE")
BMIN=10; N=15; V=43;QL=10**15;DF=10**36
A=10306;C=21954;D=100000;RU=F('0.31814293085928')

def edges(B,faces):
 assert len(faces)==28-B
 return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+15+k,B+c) for k,fa in enumerate(faces) for c in fa]

def signature(B,faces):
 adj=[1<<i for i in range(B)]
 for fa in faces:
  bs=[v for v in fa if v<B]
  for u in bs:
   for v in bs:adj[u]|=1<<v
 masks=set()
 for i in range(B):
  reach=adj[i]|adj[(i+1)%B]
  for j in range(i+1,B):
   if not(reach&((1<<j)|(1<<((j+1)%B)))):continue
   k=j-i
   if min(k,B-k)<=2:continue
   fl=k*A<D-C;fr=(B-k)*A<D-C
   if fl and fr:return None
   if fl:mask=((1<<j)-1)^((1<<i)-1)
   elif fr:mask=((1<<B)-1)^(((1<<j)-1)^((1<<i)-1))
   else:continue
   masks.add(mask)
 return tuple(sorted(masks))

def bounds(B,faces):
 M=[[10**8]*B for _ in range(B)]
 for i in range(B):M[i][i]=0
 for i in range(B-1):M[i][i+1]=A;M[i+1][i]=0
 M[B-1][0]=A-D;M[0][B-1]=D
 key=signature(B,faces);assert key is not None
 for mask in key:
  i=next(i for i in range(B) if mask>>i&1 and not(mask>>((i-1)%B)&1));j=next(i for i in range(B) if not(mask>>i&1) and mask>>((i-1)%B)&1)
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
 return nl[1:],nh[1:]
def ceilroot(n):
 r=isqrt(n);return r+(r*r<n)

def check(B,ed,lo,hi,c):
 lam=c['lambda_num'];dv=c['support_num'];flow=c['force_num'];assert len(lam)==len(dv)==B-1 and len(flow)==len(ed)
 assert all(type(x)is int and x>=0 for x in lam) and all(len(z)==2 and all(type(v)is int for v in z) for z in dv+flow)
 aa=[];bb=[]
 for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
 target=[[0,0] for _ in range(V)];act=[[0,0] for _ in range(V)]
 for i in range(B-1):
  for d in range(2):
   target[i+1][d]=(lam[i]*aa[i][d]-QL*dv[i][d])*(DF//(QL*QA));target[0][d]-=target[i+1][d]
 for (u,v),z in zip(ed,flow):
  for d in range(2):act[u][d]+=z[d];act[v][d]-=z[d]
 assert act==target,'integer force balance'
 ns=sum(ceilroot(x*x+y*y) for x,y in flow)
 bn=sum(lam[i]*(bb[i]-aa[i][0]*(QB//QA))-QL*(ceilroot(dv[i][0]**2+dv[i][1]**2)-dv[i][0])*(QB//QA) for i in range(B-1))
 if ns<=0:return None
 return F(bn*DF,QL*QB*ns)-RU
