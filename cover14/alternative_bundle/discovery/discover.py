#!/usr/bin/env python3
"""Numerical discovery only; not a proof. Uses variable circle anchors and face witnesses."""
import numpy as np, json, time
from scipy.optimize import minimize
from scipy.spatial import Delaunay
from pathlib import Path
BASE=Path(__file__).resolve().parent
N=14

def model(B,faces):
 nf=len(faces);na=B-1;nt=na+2*(N+nf);nv=nt+1
 edges=[(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+k,B+c) for k,f in enumerate(faces) for c in f]
 eu=np.array([e[0] for e in edges]);ev=np.array([e[1] for e in edges]);ne=len(edges)
 def xyz(x):
  th=np.r_[0,x[:na]];return np.vstack([np.column_stack([np.cos(th),np.sin(th)]),x[na:nt].reshape(-1,2)])
 def con(x):
  z=xyz(x);dz=z[eu]-z[ev];return x[-1]-(dz*dz).sum(1)
 def jac(x):
  z=xyz(x);dz=z[eu]-z[ev];J=np.zeros((ne,nv));J[:,-1]=1
  for nodes,sgn in [(eu,-2),(ev,2)]:
   ii=np.where(nodes>=B)[0];nn=nodes[ii];cols=na+2*(nn-B)
   J[ii,cols]+=sgn*dz[ii,0];J[ii,cols+1]+=sgn*dz[ii,1]
   ii=np.where((nodes>0)&(nodes<B))[0];nn=nodes[ii]
   J[ii,nn-1]+=sgn*(-dz[ii,0]*z[nn,1]+dz[ii,1]*z[nn,0])
  return J
 def objective(x):return x[-1]
 og=np.zeros(nv);og[-1]=1
 return xyz,con,jac,objective,lambda x:og,nt

def solve(B,C,faces=None,th=None):
 faces=np.array(sorted(tuple(sorted(map(int,f))) for f in Delaunay(C).simplices)) if faces is None else np.array(faces)
 if len(faces)!=2*N-B-2:return None
 nf=len(faces);xyz,con,jac,fun,grad,nt=model(B,faces)
 if th is None:th=2*np.pi*np.arange(1,B)/B
 ps=[]
 for f in faces:
  pts=C[f];a=2*(pts[1:]-pts[0]);b=(pts[1:]**2).sum(1)-(pts[0]**2).sum()
  try:p=np.linalg.solve(a,b)
  except:p=pts.mean(0)
  ps.append(p)
 x=np.r_[th,C.ravel(),np.array(ps).ravel(),0.13]
 x[-1]+=max(0,-con(x).min())
 bds=[(2*np.pi*k/B-.31,2*np.pi*k/B+.31) for k in range(1,B)]+[(-1.2,1.2)]*(2*(N+nf))+[(.07,.3)]
 opt=minimize(fun,x,jac=grad,method='SLSQP',constraints={'type':'ineq','fun':con,'jac':jac},bounds=bds,options={'ftol':2e-14,'maxiter':1500})
 z=xyz(opt.x)
 out={'N':N,'B':B,'faces':faces.tolist(),'centers':z[B:B+N].tolist(),'anchors':z[:B].tolist(),'points':z[B+N:].tolist(),'angles':np.r_[0,opt.x[:B-1]].tolist(),'t':opt.x[-1],'r':float(np.sqrt(opt.x[-1])),'x':opt.x.tolist(),'success':bool(opt.success),'message':str(opt.message),'min_constraint':float(con(opt.x).min()),'nit':opt.nit}
 return out

def main():
 rng=np.random.default_rng(513412);best=1.;st=time.time()
 for k in range(35):
  B=10 if k<25 else 11;I=N-B
  th=2*np.pi*(np.arange(B)-.5)/B
  ang=2*np.pi*np.arange(I)/I+(k%8)*np.pi/8
  C=np.vstack([.81*np.column_stack([np.cos(th),np.sin(th)]),.33*np.column_stack([np.cos(ang),np.sin(ang)])])
  if k%8>=4:C+=rng.normal(0,.025,C.shape)
  try:o=solve(B,C)
  except Exception as e:print('error',k,repr(e),flush=True);continue
  if o is None:continue
  print(k, o['r'],o['success'],o['min_constraint'],o['nit'],round(time.time()-st,2),flush=True)
  (BASE/f'discovery_{k}.json').write_text(json.dumps(o,indent=2))
  if o['r']<best and o['min_constraint']>-1e-9:
   best=o['r'];(BASE/'candidate_numeric.json').write_text(json.dumps(o,indent=2))
  if best<.3318:break
if __name__=='__main__':main()
