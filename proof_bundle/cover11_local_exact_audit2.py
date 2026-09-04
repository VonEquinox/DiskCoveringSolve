#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import combinations
import json,math,sys,time
sys.path.insert(0,'/mnt/data')
from cover11_fixed_trig import Q as TQ,enc,sin_iv,cos_iv
L=json.load(open('/mnt/data/cover11_local_regularized_certificate.json'));Odata=json.load(open('/mnt/data/cover11_rational_orthogonal.json'))
selected=list(map(int,L['selected']));pairs=[tuple(map(int,p)) for p in L['pairs']];x=[F(z) for z in L['gap_center']]
O=[[F(int(a),int(b)) for a,b in row] for row in Odata['matrix']];A=[[O[i][j] for j in range(4)] for i in range(8)];N=[[O[i][4+j] for j in range(4)] for i in range(8)]
# verify exact orthogonality
for i in range(8):
 for j in range(8):assert sum((O[k][i]*O[k][j] for k in range(8)),F(0))==F(i==j)
# Recompute every stored proposal Kron row from the claimed original graph
# weights.  Anchors are nodes 0,...,8; all other nodes are eliminated.  The
# factor 2 matches E(g)=sum D_ij(1-cos(theta_ij)).
edges=[tuple(map(int,e)) for e in L['edges']];weight_den=int(L['weight_den'])
assert len(edges)==45 and len(L['weight_nums'])==len(selected)
assert pairs==[(i,j) for i in range(9) for j in range(i+1,9)]
def exact_kron_D(nums):
 assert len(nums)==len(edges) and min(nums)>=0 and sum(nums)==weight_den
 nnode=1+max(max(e) for e in edges);adj={i:{} for i in range(nnode)}
 def add(a,b,z):
  if a==b or z==0:return
  adj[a][b]=adj[a].get(b,F(0))+z
  adj[b][a]=adj[b].get(a,F(0))+z
 for (a,b),num in zip(edges,nums):add(a,b,F(num,weight_den))
 free=set(range(9,nnode))
 while free:
  v=min(free,key=lambda z:(sum(1 for x in adj.get(z,{}).values() if x),z));free.remove(v)
  nei=[(u,z) for u,z in adj[v].items() if z];mass=sum((z for _,z in nei),F(0))
  assert mass>0
  for i,(a,xv) in enumerate(nei):
   for b,yv in nei[i+1:]:add(a,b,xv*yv/mass)
  for a,_ in nei:adj[a].pop(v,None)
  del adj[v]
 return [2*adj.get(i,{}).get(j,F(0)) for i,j in pairs]
for idx,nums0 in zip(selected,L['weight_nums']):
 recomputed=exact_kron_D(list(map(int,nums0)))
 stored=[F(int(a),int(b)) for a,b in L['D'][str(idx)]]
 assert recomputed==stored,(idx,'stored Kron row does not match weight_nums')
print('SIX PROPOSAL KRON ROWS RECOMPUTED EXACTLY',flush=True)

# round exact proposal conductances to central dyadics; record approximation only for diagnostics
DBITS=110;DQ=1<<DBITS;D=[];derr=F(0)
for idx in selected:
 row=[]
 for p,q in L['D'][str(idx)]:
  z=F(int(p),int(q));v=z*DQ;n=(2*v.numerator+v.denominator)//(2*v.denominator);r=F(n,DQ);derr=max(derr,abs(r-z));row.append(r)
 D.append(row)
print('D rounding error',float(derr),flush=True)
# rational central derivatives
G=[];Hs=[]
for row in D:
 g=[F(0)]*8;H=[[F(0)]*8 for _ in range(8)]
 for co,(i,j) in zip(row,pairs):
  z=sum(x[i:j],F(0));sv=sin_iv(enc(z));cv=cos_iv(enc(z));sm=F(sv.lo+sv.hi,2*TQ);cm=F(cv.lo+cv.hi,2*TQ)
  for a in range(i,j):
   g[a]+=co*sm
   for b in range(i,j):H[a][b]+=co*cm
 G.append(g);Hs.append(H)
def mt(A):return [list(x) for x in zip(*A)]
def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0)) for j in range(len(B[0]))] for i in range(len(A))]
def mv(A,v):return [sum((a*x for a,x in zip(r,v)),F(0)) for r in A]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def solve(M,b):
 n=len(M);A=[list(M[i])+[b[i]] for i in range(n)]
 for k in range(n):
  p=next(i for i in range(k,n) if A[i][k]);A[k],A[p]=A[p],A[k];z=A[k][k];A[k]=[q/z for q in A[k]]
  for i in range(n):
   if i==k:continue
   z=A[i][k]
   if z:A[i]=[A[i][j]-z*A[k][j] for j in range(n+1)]
 return [A[i][-1] for i in range(n)]
def det_nonzero(M):
 try:solve(M,[F(0)]*len(M));return True
 except StopIteration:return False
def ldlt_pd(M):
 n=len(M);L0=[[F(0)]*n for _ in range(n)];d=[F(0)]*n
 for i in range(n):
  z=M[i][i]-sum((L0[i][k]*L0[i][k]*d[k] for k in range(i)),F(0));assert z>0;d[i]=z;L0[i][i]=1
  for j in range(i+1,n):L0[j][i]=(M[j][i]-sum((L0[j][k]*L0[i][k]*d[k] for k in range(i)),F(0)))/d[i]
 return d
# Project gradients into exact rational orthonormal A,N
P=[];res=[]
At=mt(A);Nt=mt(N)
for g in G:
 P.append(mv(At,g));res.append(mv(Nt,g))
resnorm=max(sum((z*z for z in r),F(0)) for r in res)
assert resnorm<F(1,10**22) # each residual norm <1e-11
print('projection residual norm',math.sqrt(float(resnorm)),flush=True)
# sigma_min(P)>=0.007
Pt=mt(P);Gram=mm(Pt,P);s0=F(7,1000);ldlt_pd([[Gram[i][j]-(s0*s0 if i==j else 0) for j in range(4)] for i in range(4)])
# exact positive convex combination to prove polar bounded
lam6=F(763854,10**7)
Mlam=[[P[j][i] for j in range(5)] for i in range(4)]+[[F(1)]*5]
blam=[-lam6*P[5][i] for i in range(4)]+[1-lam6]
lam=solve(Mlam,blam)+[lam6];assert min(lam)>0
for k in range(4):assert sum((lam[i]*P[i][k] for i in range(6)),F(0))==0
assert sum(lam,F(0))==1
print('lambda min',float(min(lam)),flush=True)
# polar vertices and inradius >= 0.00135
kap=F(27,20000);feas=[]
for Iset in combinations(range(6),4):
 M=[P[i] for i in Iset]
 try:y=solve(M,[F(1)]*4)
 except StopIteration:continue
 vals=[dot(P[i],y) for i in range(6)]
 if max(vals)<=1:
  n2=dot(y,y);assert n2<1/(kap*kap),(Iset,float(n2),float(1/(kap*kap)))
  feas.append((Iset,n2))
assert len(feas)==8,len(feas)
print('polar vertices',len(feas),'maxnorm',math.sqrt(float(max(z for _,z in feas))),'inradius lower',float(kap),flush=True)
# Hessian block constants
mhat=F(27,10000);bhat=F(36,625);phat=F(1137,2000)
for ii,H in enumerate(Hs):
 HA=mm(mt(A),mm(H,A));HN=mm(mt(N),mm(H,N));X=mm(mt(A),mm(H,N))
 ldlt_pd([[HN[i][j]-(mhat if i==j else 0) for j in range(4)] for i in range(4)])
 XtX=mm(mt(X),X);ldlt_pd([[(bhat*bhat if i==j else F(0))-XtX[i][j] for j in range(4)] for i in range(4)])
 ldlt_pd([[(phat if i==j else F(0))-HA[i][j] for j in range(4)] for i in range(4)])
 ldlt_pd([[(phat if i==j else F(0))+HA[i][j] for j in range(4)] for i in range(4)])
 print('H block',selected[ii],'passed',flush=True)
# third derivative and full Hessian crude coefficient bounds
sqrtU={}
for n in range(1,9):
 q=10**15;p=math.isqrt(n*q*q)
 if p*p<n*q*q:p+=1
 assert p*p>=n*q*q;sqrtU[n]=F(p,q)
M3=[];Habs=[]
for row in D:
 M3.append(sum((co*(j-i)*sqrtU[j-i] for co,(i,j) in zip(row,pairs)),F(0)))
 Habs.append(sum((co*(j-i) for co,(i,j) in zip(row,pairs)),F(0)))
print('raw M3max',float(max(M3)),'raw Habsmax',float(max(Habs)),flush=True)
assert max(M3)<F(307,100);assert max(Habs)<F(3,2)
print('M3max',float(max(M3)),'Hnorm crude',float(max(Habs)),flush=True)
# final elementary radius conditions
k=F(13,10000);m=F(13,5000);b=F(29,500);p=F(57,100);M=F(307,100);rho=F(11,5000)
assert m/2-M*rho/6>0
assert k-(b+p/2)*rho-M*rho*rho/6>0
print('FINAL local inequalities',float(m/2-M*rho/6),float(k-(b+p/2)*rho-M*rho*rho/6),flush=True)
print('LOCAL CENTRAL CERTIFICATE VERIFIED')
