"""Untrusted convex optimization discovers exact integer force certificates."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,time,math
import numpy as np
from scipy.optimize import minimize,linprog
from exact_arcs import QA,QB,arc_halfplane
R=Path(__file__).resolve().parent
BMIN=12; N=19; V=55;QL=10**15;DF=10**36
A=894563;C=1871672;C6=3128331;D=10000000;RU=F('0.27735009811262')

def edges(B,faces):
 assert len(faces)==2*N-B-2
 return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+k,B+c) for k,fa in enumerate(faces) for c in fa]

def signature(B,faces):
 adj=[1<<i for i in range(N)]
 for fa in faces:
  mask=sum(1<<i for i in fa)
  for v in fa:adj[v]|=mask
 reach2=[]
 for v in range(N):
  r=adj[v]
  for u in range(N):
   if adj[v]>>u&1:r|=adj[u]
  reach2.append(r)
 masks=set()
 for i in range(B):
  r1=adj[i]|adj[(i+1)%B];r2=reach2[i]|reach2[(i+1)%B]
  for j in range(i+1,B):
   ends=(1<<j)|(1<<((j+1)%B))
   if r1&ends:cap=C;tag=0
   elif r2&ends:cap=C6;tag=1<<B
   else:continue
   k=j-i
   if min(k,B-k)*A<=cap:continue
   fl=k*A<D-cap;fr=(B-k)*A<D-cap
   if fl and fr:return None
   if fl:mask=((1<<j)-1)^((1<<i)-1)
   elif fr:mask=((1<<B)-1)^(((1<<j)-1)^((1<<i)-1))
   else:continue
   masks.add(tag|mask)
 return tuple(sorted(masks))

def bounds(B,faces):
 M=[[10**8]*B for _ in range(B)]
 for i in range(B):M[i][i]=0
 for i in range(B-1):M[i][i+1]=A;M[i+1][i]=0
 M[B-1][0]=A-D;M[0][B-1]=D
 key=signature(B,faces);assert key is not None
 for encoded in key:
  mask=encoded&((1<<B)-1);cap=C6 if encoded>>B else C
  i=next(i for i in range(B) if mask>>i&1 and not(mask>>((i-1)%B)&1));j=next(i for i in range(B) if not(mask>>i&1) and mask>>((i-1)%B)&1)
  M[i][j]=min(M[i][j],cap-(D if j<i else 0))
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

def rationalize(B,ed,lo,hi,lam,support,flow):
 la=[max(0,int(round(float(v)*QL))) for v in lam];dv=[[int(round(float(v)*QA)) for v in row] for row in support];ff=[[int(round(float(v)*DF)) for v in row] for row in flow]
 aa=[arc_halfplane(l,h)[0] for l,h in zip(lo,hi)]
 residual=[[0,0] for _ in range(V)]
 for i in range(B-1):
  for d in range(2):residual[i+1][d]=(la[i]*aa[i][d]-QL*dv[i][d])*(DF//(QL*QA));residual[0][d]-=residual[i+1][d]
 adj=[[] for _ in range(V)]
 for k,((u,v),z) in enumerate(zip(ed,ff)):
  adj[u].append((v,k));adj[v].append((u,k))
  for d in range(2):residual[u][d]-=z[d];residual[v][d]+=z[d]
 order=[0];par={0:(-1,-1)}
 for u in order:
  for v,k in adj[u]:
   if v not in par:par[v]=(u,k);order.append(v)
 assert len(order)==V
 for v in order[:0:-1]:
  u,k=par[v];sgn=1 if ed[k][0]==v else -1
  for d in range(2):ff[k][d]+=sgn*residual[v][d];residual[u][d]+=residual[v][d];residual[v][d]=0
 assert residual[0]==[0,0]
 out={'lambda_num':la,'support_num':dv,'force_num':ff}
 margin=check(B,ed,lo,hi,out)
 if margin is not None and margin>0:return out,margin
 return None,margin

class Solver:
 def __init__(self,B,faces):
  self.B=B;self.ed=edges(B,faces);self.E=np.array(self.ed);self.ne=len(self.ed)
  self.inc=np.zeros((self.ne,V-1))
  L=np.zeros((V,V))
  for k,(u,v) in enumerate(self.ed):
   if u:self.inc[k,u-1]=1
   if v:self.inc[k,v-1]=-1
   L[u,u]+=1;L[v,v]+=1;L[u,v]-=1;L[v,u]-=1
  self.harm=-np.linalg.solve(L[B:,B:],L[B:,:B]);self.nvar=2*(V-1)+1
 def solve_slsqp(self,lo,hi,maxiter=100):
  B=self.B;E=self.E;ne=self.ne;nvar=self.nvar
  aa=[];bb=[]
  for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
  aa=np.array(aa,float)/QA;bb=np.array(bb,float)/QB
  def values(x):
   z=np.vstack([[1.,0.],x[:-1].reshape(-1,2)]);ds=z[E[:,0]]-z[E[:,1]];lens=np.linalg.norm(ds,axis=1);an=np.linalg.norm(z[1:B],axis=1)
   return z,ds,lens,an
  def fun(x):
   z,ds,lens,an=values(x)
   return np.r_[x[-1]-lens,np.einsum('ij,ij->i',aa,z[1:B])-bb,1-an]
  def jac(x):
   z,ds,lens,an=values(x);direct=ds/np.maximum(lens[:,None],1e-15);J=np.zeros((ne+2*(B-1),nvar));J[:ne,-1]=1
   J[:ne,:-1]=-(self.inc[:,:,None]*direct[:,None,:]).reshape(ne,-1)
   for i in range(B-1):J[ne+i,2*i:2*i+2]=aa[i];J[ne+B-1+i,2*i:2*i+2]=-z[i+1]/max(an[i],1e-15)
   return J
  th=np.r_[0.,[2*math.pi*float((l+h)/2) for l,h in zip(lo,hi)]];Q=np.c_[np.cos(th),np.sin(th)];Z=np.vstack([Q,self.harm@Q]);x0=np.r_[Z[1:].ravel(),max(np.linalg.norm(Z[E[:,0]]-Z[E[:,1]],axis=1))+.01]
  obj=np.zeros(nvar);obj[-1]=1
  rr=minimize(lambda x:x[-1],x0,jac=lambda x:obj,constraints={'type':'ineq','fun':fun,'jac':jac},method='SLSQP',options={'ftol':1e-11,'maxiter':maxiter,'disp':False})
  z,ds,lens,an=values(rr.x);status={'success':bool(rr.success),'radius':float(rr.fun),'iterations':int(rr.nit),'violation':float(min(fun(rr.x)))}
  if hasattr(rr,'multipliers'):
   mul=rr.multipliers;lam=mul[ne:ne+B-1];support=mul[ne+B-1:,None]*z[1:B]/np.maximum(an[:,None],1e-15);flow=mul[:ne,None]*ds/np.maximum(lens[:,None],1e-15)
   cert,margin=rationalize(B,self.ed,lo,hi,lam,support,flow)
   status['exact_margin']=None if margin is None else float(margin)
   return cert,rr.x,status
  return None,rr.x,status

if __name__=='__main__':
 rr=json.load(open(R/'metric19_residuals.json'));st=time.time()
 for k in [0,1,100,931,932,1000,4000,len(rr)-1]:
  row=rr[k];B=row['B'];lo,hi=rootbox(bounds(B,row['faces']));ss=Solver(B,row['faces']);c,x,status=ss.solve_slsqp(lo,hi);print(k,B,status,'cert',c is not None,'elapsed',time.time()-st,flush=True)
