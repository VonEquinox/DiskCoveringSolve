#!/usr/bin/env python3
"""Exact anchor-only isolation, independent of right inverses and LICQ.
Every accepting inequality is checked with integers/Fraction; the congruence
is an untrusted, precomputed rational certificate.
"""
from fractions import Fraction as F
from pathlib import Path
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
import json,hashlib,sys,time
ROOT=Path(__file__).resolve().parent
import system14 as S
Q=10**15; QR=10**12; DEN=Q**4
GAMMA=300; LAMBDA=F(73,20000); NU=F(147,10000)
BUFFER=F(1,10**9); EPS=F(2,Q); RADIUS=F(1,6)

def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
 c=json.loads((ROOT/'root14.json').read_text())
 qx=int(c['Qx']); rho=F(1,10**80)
 x=[F(int(v),qx) for v in c['xnum']]
 def rnd(v):
  y=v*Q;return (2*y.numerator+y.denominator)//(2*y.denominator)
 a=[rnd(v) for v in x]
 for i in list(range(62))+list(range(63,116)):
  assert F(a[i],Q)-EPS<=x[i]-rho<=x[i]+rho<=F(a[i],Q)+EPS
 assert all(F(0)<x[i]-rho and x[i]+rho<F(1,10) for i in range(63,107))
 assert all(-NU<x[i]-rho and x[i]+rho<0 for i in range(107,116))
 coord=[[Q,0]]+[a[2*i:2*i+2] for i in range(31)]
 edges=S.EDGES
 assert len(edges)==44 and len(set(edges))==44
 A=[[0]*62 for _ in range(53)];M=[[0]*62 for _ in range(62)];degree=[0]*32
 for k,(u,v) in enumerate(edges):
  degree[u]+=1;degree[v]+=1
  for d in range(2):
   dif=coord[u][d]-coord[v][d]
   for n,sgn in ((u,1),(v,-1)):
    if n:A[k][2*(n-1)+d]+=2*sgn*dif
   for n,s in ((u,1),(v,-1)):
    if n:
     for m,t in ((u,1),(v,-1)):
      if m:M[2*(n-1)+d][2*(m-1)+d]+=a[63+k]*s*t
 for i in range(9):
  for d in range(2):
   A[44+i][2*i+d]=2*coord[i+1][d]
   M[2*i+d][2*i+d]+=a[107+i]
 # diagonal augmentation weights: w_e^2 on rods, NU^2 on circle equations.
 nuq=NU*Q;assert nuq.denominator==1
 di=[a[63+i]**2 for i in range(44)]+[nuq.numerator**2]*9
 H=[[M[i][j]*Q**3+GAMMA*sum(di[k]*A[k][i]*A[k][j] for k in range(53)) for j in range(62)] for i in range(62)]
 for i in range(62):
  sub=(BUFFER+(LAMBDA if i<18 else 0))*DEN;assert sub.denominator==1
  H[i][i]-=sub.numerator
 # Rigorous operator perturbation budget, exact root -> rational rounded matrix.
 assert sum(v*v for row in A for v in row)<81*Q*Q
 dA=236*EPS; dM=(2*max(degree)+1)*EPS
 dD=EPS/4 # |w^2-wbar^2| <= (|w|+|wbar|) EPS < EPS/4.
 assert 100*max(di)<Q*Q
 error=dM+GAMMA*(F(1,100)*(18+dA)*dA+dD*(9+dA)**2)
 assert error<BUFFER
 # Feasibility gives ||W A_rod d||_1 <= NU ||d_anchor||^2.
 # Consequently d^T A^T diag(w^2,NU^2) A d <=2 NU^2||d_anchor||^4.
 assert 2*GAMMA*NU**2*RADIUS**2<LAMBDA
 return H,error,c

def main():
 t0=time.time();H,error,c=prepare();path=ROOT/'anchor_isolation_certificate.json'
 cert=json.loads(path.read_text());assert cert['QR']==QR
 assert cert['root_sha256']==h(ROOT/'root14.json')
 assert cert['root_data_sha256']==h(ROOT/'system14.py')
 R=cert['Rnum'];assert len(R)==62 and all(len(row)==62 for row in R)
 assert all(isinstance(v,int) for row in R for v in row)
 assert all(R[i][j]==0 for i in range(62) for j in range(i))
 assert all(R[i][i]!=0 for i in range(62))
 HR=[[sum(H[i][k]*R[k][j] for k in range(j+1)) for j in range(62)] for i in range(62)]
 C=[[sum(R[k][i]*HR[k][j] for k in range(i+1)) for j in range(62)] for i in range(62)]
 assert all(C[i][j]==C[j][i] for i in range(62) for j in range(62))
 margins=[C[i][i]-sum(abs(C[i][j]) for j in range(62) if i!=j) for i in range(62)]
 assert min(margins)>0
 out={'verified':True,'theorem':'Any realization of the 44-rod candidate graph with all rods <= r_star, unit anchors, fixed q0=(1,0), and ||Q-Q_star||_F <=1/6, is the candidate itself.',
 'gamma':GAMMA,'lambda':str(LAMBDA),'nu_upper':str(NU),'anchor_Frobenius_radius':str(RADIUS),
 'per_anchor_angle_halfwidth':str(F(1,18)),
 'matrix_perturbation_bound':str(error),'matrix_buffer':str(BUFFER),
 'isolation_margin':str(LAMBDA-2*GAMMA*NU**2*RADIUS**2),
 'congruence_diagdom_margin_normalized':str(F(min(margins),DEN*QR*QR)),
 'root_sha256':cert['root_sha256'],'certificate_sha256':h(path),'seconds':time.time()-t0}
 (ROOT/'anchor_isolation_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
