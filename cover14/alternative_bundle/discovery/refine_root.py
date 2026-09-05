"""Discover a high precision full (not symmetry-restricted) KKT root and inverse."""
import numpy as np, scipy.linalg as la
from scipy.optimize import root,nnls
import mpmath as mp,json,time
from pathlib import Path
from framework import *
BASE=Path(__file__).resolve().parent
c=json.loads((BASE/'candidate_numeric.json').read_text());Z=np.vstack([c['anchors'],c['centers'],np.array(c['points'])[ACTIVE_IDS]])
x=np.r_[Z[1:].ravel(),c['t']]
z=np.r_[x,np.zeros(NC)]
f,j=values_jac(z);G=np.array(j)[:NC,:ND+1]
# Constrained rod multipliers and signed circle multipliers by least squares.
from scipy.optimize import lsq_linear
ls=lsq_linear(G.T,-np.eye(ND+1)[-1],bounds=(np.r_[np.zeros(NE),np.full(B-1,-np.inf)],np.full(NC,np.inf)),tol=1e-15,max_iter=1000)
z=np.r_[x,ls.x]
def fun(z):f,j=values_jac(z);return np.array(f,dtype=float)
def jac(z):f,j=values_jac(z);return np.array(j,dtype=float)
r=root(fun,z,jac=jac,tol=1e-11);z=r.x
print('floatroot',r.success,'res',max(abs(fun(z))),'wmin',min(z[WSTART:MUSTART]),'t',z[TID],flush=True)
f,J=values_jac(z);J=np.array(J);print('singvals',np.linalg.svd(J,compute_uv=False)[-8:],flush=True)
A=J[:NC,:ND];M=J[NC:NC+ND,:ND]/2
V=la.null_space(A);print('tangent inertia',np.linalg.eigvalsh(V.T@M@V),flush=True)
mp.mp.dps=100;zm=mp.matrix([mp.mpf(float(a)) for a in z]);st=time.time()
for k in range(10):
 f,j=values_jac(list(zm));fm=mp.matrix(f);err=mp.norm(fm,mp.inf);print('Newton',k,mp.nstr(err,8),round(time.time()-st,2),flush=True)
 if err<mp.mpf('1e-95'):break
 dz=mp.lu_solve(mp.matrix(j),-fm);zm+=dz
f,j=values_jac(list(zm));Y=mp.inverse(mp.matrix(j));print('inverse',round(time.time()-st,2),flush=True)
Qx=10**90;Qy=10**75
rnd=lambda v,Q:int(mp.nint(v*Q))
cert={'schema':'r14-kkt-rational-root-v1','N':N,'Qx':str(Qx),'Qy':str(Qy),'rho_num':'1','rho_den':str(10**65),'xnum':[str(rnd(v,Qx)) for v in zm],'Ynum':[[str(rnd(Y[i,j],Qy)) for j in range(N)] for i in range(N)]}
(BASE/'root_certificate.json').write_text(json.dumps(cert,separators=(',',':')))
rr=mp.sqrt(zm[TID]);report={'dimension':N,'r':mp.nstr(rr,95),'R':mp.nstr(1/rr,95),'t':mp.nstr(zm[TID],95),'w_min':mp.nstr(min(zm[WSTART:MUSTART]),30),'z':[mp.nstr(v,100) for v in zm]}
(BASE/'root_highprec.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='z'},indent=2),flush=True)
