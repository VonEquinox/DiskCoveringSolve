#!/usr/bin/env python3
"""Exact, anchor-only isolation theorem for the 44-rod candidate framework."""
from fractions import Fraction as F
from pathlib import Path
import json,time,hashlib
from framework import *
from verify_root import load
if not __debug__:raise RuntimeError('Exact verification requires assertions')
BASE=Path(__file__).resolve().parent
Q=10**15;QR=10**12;DEN=Q**4
GAMMA=100;LAMBDA=F(7,5000);NU=F(147,10000);RADIUS=F(7,40)
BUFFER=F(1,10**9);EPS=F(2,Q)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def prepare():
 c,x,Y,rho=load()
 def rnd(v):
  y=v*Q;return (2*y.numerator+y.denominator)//(2*y.denominator)
 a=[rnd(v) for v in x]
 for i in list(range(ND))+list(range(WSTART,N)):
  assert F(a[i],Q)-EPS<=x[i]-rho<=x[i]+rho<=F(a[i],Q)+EPS
 assert all(0<x[i]-rho and x[i]+rho<F(1,10) for i in range(WSTART,MUSTART))
 assert all(-NU<x[i]-rho and x[i]+rho<0 for i in range(MUSTART,N))
 pos=[[Q,0]]+[a[2*i:2*i+2] for i in range(NN-1)]
 A=[[0]*ND for _ in range(NC)];M=[[0]*ND for _ in range(ND)];degree=[0]*NN
 for k,(u,v) in enumerate(EDGES):
  degree[u]+=1;degree[v]+=1
  for d in range(2):
   dif=pos[u][d]-pos[v][d]
   for n,sgn in ((u,1),(v,-1)):
    if n:A[k][2*(n-1)+d]+=2*sgn*dif
   for n,s in ((u,1),(v,-1)):
    if n:
     for m,t in ((u,1),(v,-1)):
      if m:M[2*(n-1)+d][2*(m-1)+d]+=a[WSTART+k]*s*t
 for i in range(B-1):
  for d in range(2):
   A[NE+i][2*i+d]=2*pos[i+1][d]
   M[2*i+d][2*i+d]+=a[MUSTART+i]
 nuq=NU*Q;assert nuq.denominator==1
 di=[a[WSTART+i]**2 for i in range(NE)]+[nuq.numerator**2]*(B-1)
 H=[[M[i][j]*Q**3+GAMMA*sum(di[k]*A[k][i]*A[k][j] for k in range(NC)) for j in range(ND)] for i in range(ND)]
 for i in range(ND):
  sub=(BUFFER+(LAMBDA if i<2*(B-1) else 0))*DEN;assert sub.denominator==1;H[i][i]-=sub.numerator
 # A rod row has <=4 entries, each error <=4 EPS; circle row <=2, each <=2 EPS.
 assert 64**2 > NE*4*16+(B-1)*2*4
 dA=64*EPS;dM=(2*max(degree)+1)*EPS;dD=EPS/4
 assert sum(v*v for row in A for v in row)<81*Q*Q
 assert all(abs(F(a[i],Q))+abs(x[i])+rho<F(1,4) for i in range(WSTART,MUSTART))
 assert 100*max(di)<Q*Q
 error=dM+GAMMA*(F(1,100)*(18+dA)*dA+dD*(9+dA)**2)
 assert error<BUFFER
 assert 2*GAMMA*NU**2*RADIUS**2<LAMBDA
 return H,error

def verify():
 st=time.time();H,error=prepare();cert=json.loads((BASE/'anchor_isolation_certificate.json').read_text())
 assert cert['QR']==QR and cert['root_sha256']==digest(BASE/'root_certificate.json')
 R=cert['Rnum'];assert len(R)==ND and all(len(row)==ND for row in R)
 assert all(type(v) is int for row in R for v in row)
 assert all(R[i][j]==0 for i in range(ND) for j in range(i)) and all(R[i][i]!=0 for i in range(ND))
 HR=[[sum(H[i][k]*R[k][j] for k in range(j+1)) for j in range(ND)] for i in range(ND)]
 C=[[sum(R[k][i]*HR[k][j] for k in range(i+1)) for j in range(ND)] for i in range(ND)]
 assert all(C[i][j]==C[j][i] for i in range(ND) for j in range(ND))
 margins=[C[i][i]-sum(abs(C[i][j]) for j in range(ND) if i!=j) for i in range(ND)];assert min(margins)>0
 out={'verified':True,'ND':ND,'rods':NE,'gamma':GAMMA,'lambda':str(LAMBDA),'nu':str(NU),'anchor_radius':str(RADIUS),'perturbation_bound':str(error),'matrix_buffer':str(BUFFER),'isolation_margin':str(LAMBDA-2*GAMMA*NU**2*RADIUS**2),'congruence_margin':str(F(min(margins),DEN*QR**2)),'root_sha256':cert['root_sha256'],'seconds':time.time()-st}
 (BASE/'anchor_isolation_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':verify()
