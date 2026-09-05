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
import system18 as S
Q=10**15; QR=10**12; DEN=Q**4
GAMMA=300; LAMBDA=F(1,1000); NU=F(7,625)
BUFFER=F(1,10**8); EPS=F(2,Q); RADIUS=F(11,100)

def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
 c=json.loads((ROOT/'root18.json').read_text())
 qx=int(c['Qx']); rho=F(1,10**80)
 x=[F(int(v),qx) for v in c['xnum']]
 def rnd(v):
  y=v*Q;return (2*y.numerator+y.denominator)//(2*y.denominator)
 a=[rnd(v) for v in x]
 for i in list(range(S.G))+list(range(S.P,S.D)):
  assert F(a[i],Q)-EPS<=x[i]-rho<=x[i]+rho<=F(a[i],Q)+EPS
 assert all(F(0)<x[i]-rho and x[i]+rho<F(1,10) for i in range(S.P,S.P+S.M))
 assert all(-NU<x[i]-rho and x[i]+rho<0 for i in range(S.P+S.M,S.D))
 coord=[[Q,0]]+[a[2*i:2*i+2] for i in range(S.G//2)]
 edges=S.EDGES
 assert len(edges)==S.M and len(set(edges))==S.M
 A=[[0]*S.G for _ in range(S.C)];M=[[0]*S.G for _ in range(S.G)];degree=[0]*(S.G//2+1)
 for k,(u,v) in enumerate(edges):
  degree[u]+=1;degree[v]+=1
  for d in range(2):
   dif=coord[u][d]-coord[v][d]
   for n,sgn in ((u,1),(v,-1)):
    if n:A[k][2*(n-1)+d]+=2*sgn*dif
   for n,s in ((u,1),(v,-1)):
    if n:
     for m,t in ((u,1),(v,-1)):
      if m:M[2*(n-1)+d][2*(m-1)+d]+=a[S.P+k]*s*t
 for i in range(S.B-1):
  for d in range(2):
   A[S.M+i][2*i+d]=2*coord[i+1][d]
   M[2*i+d][2*i+d]+=a[S.P+S.M+i]
 # diagonal augmentation weights: w_e^2 on rods, NU^2 on circle equations.
 nuq=NU*Q;assert nuq.denominator==1
 di=[a[S.P+i]**2 for i in range(S.M)]+[nuq.numerator**2]*(S.B-1)
 H=[[M[i][j]*Q**3+GAMMA*sum(di[k]*A[k][i]*A[k][j] for k in range(S.C)) for j in range(S.G)] for i in range(S.G)]
 for i in range(S.G):
  sub=(BUFFER+(LAMBDA if i<2*(S.B-1) else 0))*DEN;assert sub.denominator==1
  H[i][i]-=sub.numerator
 # Rigorous operator perturbation budget, exact root -> rational rounded matrix.
 assert sum(v*v for row in A for v in row)<100*Q*Q
 assert S.C*S.G<10000 # At most 4 EPS error per Jacobian entry.
 dA=400*EPS; dM=(2*max(degree)+1)*EPS
 dD=EPS/4 # |w^2-wbar^2| <= (|w|+|wbar|) EPS < EPS/4.
 assert 100*max(di)<Q*Q
 error=dM+GAMMA*(F(1,100)*(20+dA)*dA+dD*(10+dA)**2)
 assert error<BUFFER
 # Feasibility gives ||W A_rod d||_1 <= NU ||d_anchor||^2.
 # Consequently d^T A^T diag(w^2,NU^2) A d <=2 NU^2||d_anchor||^4.
 assert 2*GAMMA*NU**2*RADIUS**2<LAMBDA
 return H,error,c

def main():
 t0=time.time();H,error,c=prepare();path=ROOT/'anchor_isolation_certificate.json'
 cert=json.loads(path.read_text());assert cert['QR']==QR
 assert cert['root_sha256']==h(ROOT/'root18.json')
 assert cert['root_data_sha256']==h(ROOT/'system18.py')
 R=cert['Rnum'];assert len(R)==S.G and all(len(row)==S.G for row in R)
 assert all(isinstance(v,int) for row in R for v in row)
 assert all(R[i][j]==0 for i in range(S.G) for j in range(i))
 assert all(R[i][i]!=0 for i in range(S.G))
 HR=[[sum(H[i][k]*R[k][j] for k in range(j+1)) for j in range(S.G)] for i in range(S.G)]
 C=[[sum(R[k][i]*HR[k][j] for k in range(i+1)) for j in range(S.G)] for i in range(S.G)]
 assert all(C[i][j]==C[j][i] for i in range(S.G) for j in range(S.G))
 margins=[C[i][i]-sum(abs(C[i][j]) for j in range(S.G) if i!=j) for i in range(S.G)]
 assert min(margins)>0
 out={'verified':True,'theorem':'Any realization of the 72-rod candidate graph with all rods <= r_star, unit anchors, fixed q0=(1,0), and ||Q-Q_star||_F <=11/100, is the candidate itself.',
 'gamma':GAMMA,'lambda':str(LAMBDA),'nu_upper':str(NU),'anchor_Frobenius_radius':str(RADIUS),
 'matrix_perturbation_bound':str(error),'matrix_buffer':str(BUFFER),
 'isolation_margin':str(LAMBDA-2*GAMMA*NU**2*RADIUS**2),
 'congruence_diagdom_margin_normalized':str(F(min(margins),DEN*QR*QR)),
 'root_sha256':cert['root_sha256'],'certificate_sha256':h(path),'seconds':time.time()-t0}
 (ROOT/'anchor_isolation_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
