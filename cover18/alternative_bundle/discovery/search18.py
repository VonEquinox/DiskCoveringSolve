import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import discover18_B12 as d
import numpy as np,json,sys,time
from pathlib import Path
from scipy.spatial import Delaunay
B=int(sys.argv[1]);N=18;d.B=B;d.N=N
count=int(sys.argv[2]) if len(sys.argv)>2 else 500
rng=np.random.default_rng(1880+B);path=Path(__file__).parent;best=1.;t=time.time();bestp=None
for k in range(count):
 if bestp is not None and k%3:
  P=bestp.copy();amp=[.015,.035,.06,.09,.14,.20][(k//3)%6];P[B:]+=rng.normal(0,amp,P[B:].shape)
  if k%6<3:P[:B]+=rng.normal(0,amp/3,P[:B].shape)
 else:
  th=np.arange(B)*2*np.pi/B
  th2=np.arange(N-B)*2*np.pi/(N-B)+rng.uniform(0,6.28)
  ins=.35*np.c_[np.cos(th2),np.sin(th2)]+rng.normal(0,.08,(N-B,2))
  if k%9==0:ins[0]=0;ins[1:]=.42*np.c_[np.cos(th2[1:]),np.sin(th2[1:])]
  P=np.r_[.83*np.c_[np.cos(th),np.sin(th)],ins]
 try:
  fa,E=d.data(P);r=d.solve(P,fa,E)
 except Exception as e:continue
 if r['t']<best and min(r['rod_slacks'])>-1e-9:
  best=r['t'];bestp=np.array(r['points']);(path/f'best18_B{B}.json').write_text(json.dumps(r,indent=2));print('BEST',B,k,r['r'],'R',1/r['r'],'sec',time.time()-t,flush=True)
 if k%100==0:print('PROGRESS',B,k,'sec',time.time()-t,flush=True)
print('FINISHED',B,best**.5,time.time()-t,flush=True)
