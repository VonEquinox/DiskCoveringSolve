"""Optional faster discovery. Does not change the exact acceptance kernel."""
import numpy as np,math
import force18 as base
from force18 import *
from barrier18 import solve
class Solver(base.Solver):
 def __init__(self,B,faces):
  super().__init__(B,faces)
  adj=[[] for _ in range(V)]
  for k,(u,v) in enumerate(self.ed):adj[u].append((v,k));adj[v].append((u,k))
  self.parent=np.full((V,2),-1,dtype=np.int64);order=[0]
  for u in order:
   for v,k in adj[u]:
    if v and self.parent[v,0]<0:self.parent[v]=u,k;order.append(v)
  self.order=np.array(order,dtype=np.int64);assert len(order)==V
 def solve_slsqp(self,lo,hi,maxiter=100):
  B=self.B;aa=[];bb=[]
  for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
  aa=np.array(aa,float)/QA;bb=np.array(bb,float)/QB
  if min(float(h-l) for l,h in zip(lo,hi))<1e-6:return super().solve_slsqp(lo,hi,maxiter)
  angles=np.array([2*math.pi*float((l+h)/2) for l,h in zip(lo,hi)])
  Q=np.vstack([[1.,0.],np.c_[np.cos(angles),np.sin(angles)]*((1+bb/np.linalg.norm(aa,axis=1))/2)[:,None]])
  Z=np.vstack([Q,self.harm@Q]);x0=np.r_[Z[1:].ravel(),max(np.linalg.norm(Z[self.E[:,0]]-Z[self.E[:,1]],axis=1))+.02]
  try:
   x,flow,lam,support,low,its=solve(x0,self.E,aa,bb,self.parent,self.order,float(RU))
   cert,mar=rationalize(B,self.ed,lo,hi,lam,support,flow)
   status={'success':low>float(RU) or x[-1]<float(RU),'radius':float(x[-1]),'iterations':int(its),'lower_float':low,'exact_margin':None if mar is None else float(mar),'algorithm':'logbarrier'}
   if cert is not None or x[-1]<float(RU)-1e-6:return cert,x,status
  except Exception as e:pass
  return super().solve_slsqp(lo,hi,maxiter)
