"""Untrusted numerical discovery of exact, independently checkable force certificates."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from exact_arcs import arc_halfplane,QA,QB
QL=10**15;DF=10**36;RU=F(33173203427624,10**14)
N=14;NN=40;ND=2*(NN-1)

def ceilroot(n):
 r=isqrt(n);return r+(r*r<n)

def edges(B,faces):return [(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+k,B+c) for k,fa in enumerate(faces) for c in fa]

def exact_flow(B,ed,lo,hi,mults,z):
 ne=len(ed);lam=[max(0,int(round(a*QL))) for a in mults[ne:ne+B-1]];eta=[max(0,int(round(2*a*QL))) for a in mults[ne+B-1:]]
 dv=[[int(round(float(a)*QA)) for a in zz] for zz in z[1:B]]
 aa=[];bb=[]
 for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
 target=[[0,0] for _ in range(NN)]
 for i in range(B-1):
  for d in range(2):
   target[i+1][d]=(lam[i]*aa[i][d]-eta[i]*dv[i][d])*(DF//(QL*QA));target[0][d]-=target[i+1][d]
 flow=[]
 for k,(u,v) in enumerate(ed):flow.append([int(round(2*max(0,float(mults[k]))*float(z[u,d]-z[v,d])*DF)) for d in range(2)])
 act=[[0,0] for _ in range(NN)]
 for (u,v),f in zip(ed,flow):
  for d in range(2):act[u][d]+=f[d];act[v][d]-=f[d]
 deficit=[[t-a for t,a in zip(tt,ac)] for tt,ac in zip(target,act)]
 adj=[[] for _ in range(NN)]
 for k,(u,v) in enumerate(ed):adj[u].append((v,k,1));adj[v].append((u,k,-1))
 par={0:None};order=[0]
 for u in order:
  for v,k,sgn in adj[u]:
   if v not in par:par[v]=(u,k,-sgn);order.append(v)
 assert len(order)==NN
 for v in order[:0:-1]:
  u,k,sgn=par[v]
  for d in range(2):flow[k][d]+=sgn*deficit[v][d];deficit[u][d]+=deficit[v][d];deficit[v][d]=0
 assert deficit[0]==[0,0]
 ns=sum(ceilroot(f[0]**2+f[1]**2) for f in flow)
 bn=sum(lam[i]*(bb[i]-aa[i][0]*(QB//QA))-eta[i]*(ceilroot(dv[i][0]**2+dv[i][1]**2)-dv[i][0])*(QB//QA) for i in range(B-1))
 margin=F(bn*DF,QL*QB*ns)-RU if ns else F(-1)
 rec={'kind':'DUAL','lambda_num':lam,'eta_num':eta,'support_vector_num':dv,'force_num':flow}
 return margin,rec

class Solver:
 def __init__(self,B,faces):
  self.B=B;self.faces=faces;ed=self.ed=edges(B,faces);self.eu=np.array([u for u,v in ed]);self.ev=np.array([v for u,v in ed]);self.ne=len(ed);self.nv=ND+1
  self.og=np.eye(self.nv)[-1]
  L=np.zeros((NN,NN))
  for u,v in ed:L[u,u]+=1;L[v,v]+=1;L[u,v]-=1;L[v,u]-=1
  self.harm=-np.linalg.solve(L[B:,B:],L[B:,:B]);self.last=None
 def solve(self,lo,hi,x0=None):
  B=self.B;ne=self.ne;ed=self.ed;eu=self.eu;ev=self.ev;nv=self.nv
  planes=[arc_halfplane(l,h) for l,h in zip(lo,hi)];aa=np.array([a for a,b in planes],float)/QA;bb=np.array([b for a,b in planes],float)/QB
  def xyz(x):return np.vstack([[1.,0.],x[:ND].reshape(-1,2)])
  def con(x):
   z=xyz(x);dz=z[eu]-z[ev];q=z[1:B]
   return np.r_[x[-1]-(dz*dz).sum(1),(aa*q).sum(1)-bb,1-(q*q).sum(1)]
  def jac(x):
   z=xyz(x);dz=z[eu]-z[ev];J=np.zeros((ne+2*(B-1),nv));J[:ne,-1]=1
   for nodes,sgn in [(eu,-2),(ev,2)]:
    ii=np.where(nodes>0)[0];nn=nodes[ii];cols=2*(nn-1)
    J[ii,cols]+=sgn*dz[ii,0];J[ii,cols+1]+=sgn*dz[ii,1]
   rr=np.arange(B-1);J[ne+rr,2*rr]=aa[:,0];J[ne+rr,2*rr+1]=aa[:,1]
   J[ne+B-1+rr,2*rr]=-2*z[1:B,0];J[ne+B-1+rr,2*rr+1]=-2*z[1:B,1]
   return J
  if x0 is None:
   mid=np.array([float((l+h)/2)*2*np.pi for l,h in zip(lo,hi)]);q=np.vstack([[1.,0.],np.column_stack([np.cos(mid),np.sin(mid)])]);z=np.vstack([q,self.harm@q]);x0=np.r_[z[1:].ravel(),.2]
   x0[-1]=.2-min(con(x0)[:ne].min(),0.)
  out=minimize(lambda x:x[-1],x0,jac=lambda x:self.og,method='SLSQP',constraints={'type':'ineq','fun':con,'jac':jac},options={'ftol':5e-13,'maxiter':150})
  z=xyz(out.x);mult=out.multipliers
  margin,record=exact_flow(B,ed,lo,hi,mult,z)
  self.last=out.x
  return margin,record,out,z
