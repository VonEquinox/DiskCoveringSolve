import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import discover20 as d
import numpy as np,json,sys,time
from pathlib import Path
B=int(sys.argv[1]);N=20;d.B=B;d.N=N
count=int(sys.argv[2]) if len(sys.argv)>2 else 120
seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
rng=np.random.default_rng(2020+B*100+seed);path=Path(__file__).parent;best=1.;t=time.time();bestp=None
for k in range(count):
 if bestp is not None and k%3:
  P=bestp.copy();amp=[.008,.015,.028,.045,.07,.10,.14][(k//3)%7];P[B:]+=rng.normal(0,amp,P[B:].shape)
  if k%6<3:P[:B]+=rng.normal(0,amp/3,P[:B].shape)
 else:
  th=np.arange(B)*2*np.pi/B
  if k%6<3:
   th2=np.arange(N-B-1)*2*np.pi/(N-B-1)+rng.uniform(0,6.28)
   ins=np.r_[np.zeros((1,2)),.46*np.c_[np.cos(th2),np.sin(th2)]]+rng.normal(0,.025,(N-B,2))
  else:
   th2=np.arange(N-B)*2*np.pi/(N-B)+rng.uniform(0,6.28)
   ins=.37*np.c_[np.cos(th2),np.sin(th2)]+rng.normal(0,.06,(N-B,2))
  P=np.r_[.85*np.c_[np.cos(th),np.sin(th)],ins]
 try:
  fa,E=d.data(P);r=d.solve(P,fa,E)
 except Exception as e:continue
 if r['t']<best and min(r['rod_slacks'])>-1e-9:
  best=r['t'];bestp=np.array(r['points']);(path/f'best20_B{B}_s{seed}.json').write_text(json.dumps(r,indent=2));print('BEST',B,k,r['r'],'R',1/r['r'],'sec',time.time()-t,flush=True)
 if k%10==0:print('PROGRESS',B,k,'sec',time.time()-t,flush=True)
print('FINISHED',B,best**.5,time.time()-t,flush=True)
