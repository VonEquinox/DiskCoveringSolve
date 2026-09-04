#!/usr/bin/env python3
"""Generate and verify a rational norm-Krawczyk certificate for the full
110-variable active-graph KKT point of the n=11 candidate.

The system is polynomial of degree at most two.  Verification uses Fraction
arithmetic only; high-precision Newton is used solely to propose the box and
an approximate inverse.
"""
from __future__ import annotations
from fractions import Fraction as F
import json, math, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp

CERT=Path('/mnt/data/cover11_full_kkt_certificate.json')
LOG=Path('/mnt/data/cover11_full_kkt_certificate.log')
DOUBLE=Path('/mnt/data/cover11_full_kkt_root_double.json')
JOLD=json.load(open('/mnt/data/cover11_reconstructed_stresses.json'))
EDGES=[tuple(map(int,e)) for e in JOLD['edges']]
QVARS=[i for i in range(9) if i!=4]
NODEVAR={};cc=0
for i in QVARS:NODEVAR[i]=(cc,cc+1);cc+=2
for i in range(9,29):NODEVAR[i]=(cc,cc+1);cc+=2
TCOL=cc;cc+=1
NPR=cc;NCON=8+len(EDGES);NT=NPR+NCON
assert (NPR,NCON,NT)==(57,53,110)

# Constant Hessians of constraints, represented sparsely as dictionaries.
# constraints: circle rows, then edge rows
CH=[]
for node in QVARS:
 a,b=NODEVAR[node];CH.append({(a,a):2,(b,b):2})
for u,v in EDGES:
 H={}
 for d in range(2):
  iu=NODEVAR[u][d] if u in NODEVAR else None
  iv=NODEVAR[v][d] if v in NODEVAR else None
  if iu is not None:H[(iu,iu)]=H.get((iu,iu),0)+2
  if iv is not None:H[(iv,iv)]=H.get((iv,iv),0)+2
  if iu is not None and iv is not None:
   H[(iu,iv)]=H.get((iu,iv),0)-2;H[(iv,iu)]=H.get((iv,iu),0)-2
 CH.append(H)
assert len(CH)==NCON

def pos_from_y(y, zero):
 p=[[zero,zero] for _ in range(29)];p[4]=[-1+zero,zero]
 for node,(a,b) in NODEVAR.items():p[node]=[y[a],y[b]]
 return p

def con_jac(y, zero):
 p=pos_from_y(y,zero);t=y[TCOL]
 c=[zero for _ in range(NCON)];J=[[zero for _ in range(NPR)] for _ in range(NCON)];r=0
 for node in QVARS:
  a,b=NODEVAR[node];x0,x1=p[node];c[r]=x0*x0+x1*x1-1;J[r][a]=2*x0;J[r][b]=2*x1;r+=1
 for u,v in EDGES:
  d0=p[u][0]-p[v][0];d1=p[u][1]-p[v][1];c[r]=d0*d0+d1*d1-t
  if u in NODEVAR:
   a,b=NODEVAR[u];J[r][a]+=2*d0;J[r][b]+=2*d1
  if v in NODEVAR:
   a,b=NODEVAR[v];J[r][a]-=2*d0;J[r][b]-=2*d1
  J[r][TCOL]=-1;r+=1
 return c,J

def full_fun_jac(z, zero):
 y=z[:NPR];lam=z[NPR:];c,J=con_jac(y,zero)
 grad=[zero for _ in range(NPR)];grad[TCOL]=1+zero
 stat=[grad[j]+sum(J[i][j]*lam[i] for i in range(NCON)) for j in range(NPR)]
 # full Jacobian rows: constraints then stationarity
 A=[[zero for _ in range(NT)] for _ in range(NT)]
 for i in range(NCON):
  for j in range(NPR):A[i][j]=J[i][j]
 # H_L block
 for i,H in enumerate(CH):
  li=lam[i]
  for (a,b),coef in H.items():A[NCON+a][b]+=coef*li
 for j in range(NPR):
  for i in range(NCON):A[NCON+j][NPR+i]=J[i][j]
 return c+stat,A

def infnorm_mat(A):return max(sum(abs(x) for x in row) for row in A)

def qmatmul(A,B):
 nr=len(A);nk=len(A[0]);nc=len(B[0]);BT=list(zip(*B))
 return [[sum((A[i][k]*BT[j][k] for k in range(nk)),F(0)) for j in range(nc)] for i in range(nr)]

def generate():
 d=json.load(open(DOUBLE));z0=[mp.mpf(s) for s in d['z']]
 mp.mp.dps=100
 # High precision Newton.
 z=mp.matrix(z0)
 for it in range(12):
  fv,JM=full_fun_jac(list(z),mp.mpf('0'));fv=mp.matrix(fv);JM=mp.matrix(JM)
  dz=mp.lu_solve(JM,-fv);z+=dz
  md=max(abs(x) for x in dz);mr=max(abs(x) for x in fv)
  print('newton',it,'step',mp.nstr(md,8),'res',mp.nstr(mr,8),flush=True)
  if md<mp.mpf('1e-88'):break
 fv,JM=full_fun_jac(list(z),mp.mpf('0'));res=max(abs(x) for x in fv);assert res<mp.mpf('1e-80')
 print('inverse...',flush=True);Ymp=mp.inverse(mp.matrix(JM))
 # Dyadic rational midpoint and inverse.
 QB=300;Q=1<<QB
 def rnd(v):return F(int(mp.nint(v*Q)),Q)
 xq=[rnd(v) for v in z];Y=[[rnd(Ymp[i,j]) for j in range(NT)] for i in range(NT)]
 F0,J0=full_fun_jac(xq,F(0))
 print('exact residual/inverse products...',flush=True)
 shift=[sum((Y[i][j]*F0[j] for j in range(NT)),F(0)) for i in range(NT)]
 shiftinf=max(map(abs,shift));Ynorm=infnorm_mat(Y)
 YJ=qmatmul(Y,J0)
 E=[[F(i==j)-YJ[i][j] for j in range(NT)] for i in range(NT)]
 enorm=infnorm_mat(E)
 # Exact common-radius Lipschitz constant for J(z)-J(z0).
 # Compute coefficient absolute row sums directly from the quadratic structure.
 # Top constraint rows: their gradients change according to CH_i.
 Lrows=[sum(abs(v) for v in H.values()) for H in CH]
 # Bottom stationarity row a:
 #   y-block H_L changes with each multiplier: sum_{b,i}|H_i[a,b]| r
 #   lambda-block J^T changes with y: sum_{i,k}|H_i[a,k]| r
 # These two sums are equal, hence factor 2.
 for a in range(NPR):
  s=sum(abs(coef) for H in CH for (u,v),coef in H.items() if u==a)
  Lrows.append(2*s)
 L=max(Lrows);assert len(Lrows)==NT
 # common radius 1e-45 gives enormous slack
 rad=F(1,10**45)
 ratio=shiftinf/rad+enorm+Ynorm*L*rad
 print('shift/r',float(shiftinf/rad),'enorm',float(enorm),'Ynorm',float(Ynorm),'L',L,'ratio',float(ratio),flush=True)
 assert ratio<F(1,1000),ratio
 # Positivity and reduced useful facts at midpoint (root deviation is rad).
 edge_mid=xq[NPR+8:]
 edge_lower=min(edge_mid)-rad
 assert edge_lower>F(29,10000) # >0.0029
 tmid=xq[TCOL]
 cert={
  'version':1,'nprimal':NPR,'nconstraints':NCON,'total':NT,'t_index':TCOL,
  'qvars':QVARS,'nodevar':{str(k):list(v) for k,v in NODEVAR.items()},'edges':[list(e) for e in EDGES],
  'bits':QB,'midpoint':[[str(v.numerator),str(v.denominator)] for v in xq],
  'inverse':[[[str(v.numerator),str(v.denominator)] for v in row] for row in Y],
  'radius':[str(rad.numerator),str(rad.denominator)],
  'shift_inf':[str(shiftinf.numerator),str(shiftinf.denominator)],
  'inverse_residual_inf':[str(enorm.numerator),str(enorm.denominator)],
  'inverse_inf':[str(Ynorm.numerator),str(Ynorm.denominator)],
  'jacobian_lipschitz_inf':L,
  'krawczyk_ratio':[str(ratio.numerator),str(ratio.denominator)],
  'edge_multiplier_lower':[str(edge_lower.numerator),str(edge_lower.denominator)],
  't_midpoint':[str(tmid.numerator),str(tmid.denominator)],
 }
 CERT.write_text(json.dumps(cert,separators=(',',':')))
 LOG.write_text('\n'.join([
  'residual='+mp.nstr(res,30),'radius='+str(rad),'shift_over_radius='+repr(float(shiftinf/rad)),
  'inverse_residual_inf='+repr(float(enorm)),'inverse_inf='+repr(float(Ynorm)),'jacobian_lipschitz='+str(L),
  'krawczyk_ratio='+repr(float(ratio)),'edge_multiplier_lower='+repr(float(edge_lower)),
  't='+mp.nstr(z[TCOL],90),'r='+mp.nstr(mp.sqrt(z[TCOL]),90),
 ])+'\n')
 print(LOG.read_text());print(CERT,'bytes',CERT.stat().st_size)

def verify():
 c=json.load(open(CERT));assert c['total']==NT and c['nprimal']==NPR and c['nconstraints']==NCON
 xq=[F(int(a),int(b)) for a,b in c['midpoint']];Y=[[F(int(a),int(b)) for a,b in row] for row in c['inverse']]
 rad=F(int(c['radius'][0]),int(c['radius'][1]));F0,J0=full_fun_jac(xq,F(0))
 shift=[sum((Y[i][j]*F0[j] for j in range(NT)),F(0)) for i in range(NT)];shiftinf=max(map(abs,shift))
 Ynorm=infnorm_mat(Y);YJ=qmatmul(Y,J0);enorm=infnorm_mat([[F(i==j)-YJ[i][j] for j in range(NT)] for i in range(NT)])
 Lrows=[sum(abs(v) for v in H.values()) for H in CH]
 for a in range(NPR):
  s=sum(abs(coef) for H in CH for (u,v),coef in H.items() if u==a);Lrows.append(2*s)
 L=max(Lrows);ratio=shiftinf/rad+enorm+Ynorm*L*rad
 assert [str(shiftinf.numerator),str(shiftinf.denominator)]==c['shift_inf']
 assert [str(enorm.numerator),str(enorm.denominator)]==c['inverse_residual_inf']
 assert [str(Ynorm.numerator),str(Ynorm.denominator)]==c['inverse_inf']
 assert L==c['jacobian_lipschitz_inf']
 assert [str(ratio.numerator),str(ratio.denominator)]==c['krawczyk_ratio'] and ratio<F(1,1000)
 edge_lower=min(xq[NPR+8:])-rad;assert edge_lower>F(29,10000)
 print('FULL KKT RATIONAL CERTIFICATE VERIFIED')
 print('ratio',float(ratio),'edge multiplier lower',float(edge_lower),'t interval',float(xq[TCOL]-rad),float(xq[TCOL]+rad))

if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='verify':verify()
 else:generate()
