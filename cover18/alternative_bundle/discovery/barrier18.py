"""Untrusted compiled primal log-barrier SOCP discovery.
Only exact_force18.check can accept a resulting certificate.
"""
import numpy as np
from numba import njit
from math import log,sqrt
@njit(cache=True)
def objgrad(x,E,aa,bb,mu,need_h=True):
 n=len(x);r=x[-1];B=len(bb)+1;H=np.zeros((n,n));g=np.zeros(n);g[-1]=1.;f=r
 if r<=0:return np.inf,g,H
 for k in range(len(E)):
  u,v=E[k];dx=(x[2*(u-1)] if u else 1.)-(x[2*(v-1)] if v else 1.);dy=(x[2*(u-1)+1] if u else 0.)-(x[2*(v-1)+1] if v else 0.)
  s=r*r-dx*dx-dy*dy
  if s<=0:return np.inf,g,H
  f-=mu*log(s)
  if not need_h:continue
  D=np.array([dx,dy]);h1=2*mu/s;h2=4*mu/(s*s)
  g[-1]-=h1*r;H[-1,-1]+=-h1+h2*r*r
  for side in range(2):
   a=u if side==0 else v;sg=1 if side==0 else -1
   if a==0:continue
   for d in range(2):
    i=2*(a-1)+d;g[i]+=sg*h1*D[d];val=-sg*h2*r*D[d];H[i,-1]+=val;H[-1,i]+=val
    for side2 in range(2):
     b=u if side2==0 else v;tg=1 if side2==0 else -1
     if b==0:continue
     for c in range(2):
      j=2*(b-1)+c;H[i,j]+=sg*tg*(h2*D[d]*D[c]+(h1 if c==d else 0.))
 for k in range(B-1):
  qx=x[2*k];qy=x[2*k+1];v=aa[k,0]*qx+aa[k,1]*qy-bb[k];u=1-qx*qx-qy*qy
  if v<=0 or u<=0:return np.inf,g,H
  f-=mu*(log(v)+log(u))
  if not need_h:continue
  q=np.array([qx,qy]);g[2*k:2*k+2]+=2*mu*q/u-mu*aa[k]/v
  for a in range(2):
   for b in range(2):H[2*k+a,2*k+b]+=mu*aa[k,a]*aa[k,b]/(v*v)+4*mu*q[a]*q[b]/(u*u)+(2*mu/u if a==b else 0.)
 return f,g,H
@njit(cache=True)
def forces(x,E,aa,bb,mu,parent,order):
 r=x[-1];B=len(bb)+1;V=(len(x)-1)//2+1
 flow=np.zeros((len(E),2));lam=np.zeros(B-1);support=np.zeros((B-1,2));res=np.zeros((V,2))
 for k in range(B-1):
  q=x[2*k:2*k+2];lam[k]=mu/(np.dot(aa[k],q)-bb[k]);support[k]=2*mu*q/(1-np.dot(q,q));res[k+1]=lam[k]*aa[k]-support[k];res[0]-=res[k+1]
 for k in range(len(E)):
  u,v=E[k];d=np.zeros(2)
  for j in range(2):d[j]=(x[2*(u-1)+j] if u else (1. if j==0 else 0.))-(x[2*(v-1)+j] if v else (1. if j==0 else 0.))
  flow[k]=2*mu*d/(r*r-np.dot(d,d));res[u]-=flow[k];res[v]+=flow[k]
 for j in range(len(order)-1,0,-1):
  v=order[j];u,k=parent[v];sgn=1 if E[k,0]==v else -1
  flow[k]+=sgn*res[v];res[u]+=res[v];res[v]*=0
 ns=0.;bn=0.
 for k in range(len(E)):ns+=np.sqrt(np.dot(flow[k],flow[k]))
 for k in range(B-1):bn+=lam[k]*(bb[k]-aa[k,0])-np.sqrt(np.dot(support[k],support[k]))+support[k,0]
 return flow,lam,support,bn/ns
@njit(cache=True)
def solve(x,E,aa,bb,parent,order,target):
 its=0;mu=.002
 for outer in range(16):
  for it in range(55):
   its+=1;f,g,H=objgrad(x,E,aa,bb,mu)
   if not np.isfinite(f):return x,np.zeros((len(E),2)),np.zeros(len(bb)),np.zeros((len(bb),2)),-1.,its
   step=np.linalg.solve(H,-g);gd=np.dot(g,step)
   if gd>0 or not np.isfinite(gd):break
   if -gd<max(mu*1e-8,1e-17):break
   alpha=1.
   for ls in range(55):
    xx=x+alpha*step;ff,_,_=objgrad(xx,E,aa,bb,mu,False)
    if ff<=f+.01*alpha*gd:break
    alpha*=.5
   if ls==54:break
   x=xx
  flow,lam,support,lower=forces(x,E,aa,bb,mu,parent,order)
  if lower>target+5e-8 or x[-1]<target-1e-6:return x,flow,lam,support,lower,its
  mu*=.2
 return x,flow,lam,support,lower,its
