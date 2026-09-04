#!/usr/bin/env python3
"""Exact scalar transfer audit for the candidate local-minimum theorem.

Together with:
  * cover11_full_kkt_certificate.py verify
  * cover11_local_exact_audit2.py
this verifies all finite rational inequalities used by the local theorem.
The remaining steps are the explicitly stated elementary linear-algebra and
Taylor lemmas in cover11_local_theorem.md.
"""
from fractions import Fraction as F
import json,math,subprocess,sys
sys.path.insert(0,'/mnt/data')
from cover11_fixed_trig import Q as TQ,enc,sin_iv,cos_iv
import cover11_candidate_exact_core as core

# ---------- load exact root box and graph ----------
K=json.load(open('/mnt/data/cover11_full_kkt_certificate.json'))
x=[F(int(a),int(b)) for a,b in K['midpoint']]
rad=F(int(K['radius'][0]),int(K['radius'][1]))
NPR=int(K['nprimal']);TCOL=int(K['t_index']);nv={int(k):tuple(v) for k,v in K['nodevar'].items()};edges=[tuple(e) for e in K['edges']]
assert rad==F(1,10**45)
# target is strictly above exact t root box
TARGET=F(1443877286317993,10**16)
assert TARGET>x[TCOL]+rad
# all exact KKT edge multipliers positive, and stronger stored lower bound
assert min(x[NPR+8:])-rad>F(29,10000)

# ---------- equilibrium matrix B at rational midpoint ----------
pos=[]
for node in range(29):
 if node==4:pos.append((F(-1),F(0)))
 else:
  a,b=nv[node];pos.append((x[a],x[b]))
free=list(range(9,29));fm={v:k for k,v in enumerate(free)}
B0=[[F(0) for _ in range(45)] for _ in range(41)]
for e,(u,v) in enumerate(edges):
 if u in fm:
  k=fm[u];B0[2*k][e]+=pos[u][0]-pos[v][0];B0[2*k+1][e]+=pos[u][1]-pos[v][1]
 if v in fm:
  k=fm[v];B0[2*k][e]+=pos[v][0]-pos[u][0];B0[2*k+1][e]+=pos[v][1]-pos[u][1]
for e in range(45):B0[40][e]=1

def infnorm(A):return max(sum((abs(z) for z in row),F(0)) for row in A)
BT0=[list(r) for r in zip(*B0)]
assert infnorm(B0)==45 and infnorm(BT0)<3
M0=[[sum((B0[i][k]*B0[j][k] for k in range(45)),F(0)) for j in range(41)] for i in range(41)]
I=json.load(open('/mnt/data/cover11_full_B_inverse.json'));Q=1<<int(I['bits']);Y=[[F(int(I['inverse_nums'][i][j]),Q) for j in range(41)] for i in range(41)]
E=[[F(i==j)-sum((Y[i][k]*M0[k][j] for k in range(41)),F(0)) for j in range(41)] for i in range(41)]
Yinf=infnorm(Y);Einf=infnorm(E)
assert Yinf<1072 and Einf<F(1,10**12)
# Across the full KKT root box: ||dB||inf<=90 rad, ||dB^T||inf<=8 rad,
# hence ||dM||inf<700 rad.  Neumann gives ||M*^-1||inf<1100.
dM=700*rad;eta=Einf+Yinf*dM
assert eta<F(1,10**10)
Minv_bound=F(1100)
assert Yinf/(1-eta)<Minv_bound

# ---------- positive exact corrections of the six proposal stresses ----------
L=json.load(open('/mnt/data/cover11_local_regularized_certificate.json'));den=int(L['weight_den']);W0=[[F(int(n),den) for n in row] for row in L['weight_nums']]
bvec=[F(0)]*40+[F(1)]
mid_res=F(0)
for w in W0:
 assert sum(w,F(0))==1 and min(w)>F(7,100000)
 rr=[bvec[i]-sum((B0[i][j]*w[j] for j in range(45)),F(0)) for i in range(41)]
 mid_res=max(mid_res,max(map(abs,rr)))
assert mid_res<F(2,10**41)
# Exact-root residual adds at most 2*rad, since each varying edge coefficient
# changes by <=2rad and each proposal is a probability vector.
root_res=mid_res+2*rad
assert root_res<F(3,10**41)
# d w = B*^T (B*B*^T)^-1 (b-B*w0); ||B*^T||inf<4.
dwinf=4*Minv_bound*root_res
assert dwinf<F(2,10**37)
assert min(min(w) for w in W0)-dwinf>0
# Graph-energy Lipschitz lemma implies every Kron D coefficient changes by
# at most 8||dw||_1.
dD_correction=8*45*dwinf
assert dD_correction<F(1,10**34)

# ---------- exact candidate gaps are close to the rational center ----------
g=[F(z) for z in L['gap_center']]
ss=[sum(g[0:4],F(0)),sum(g[1:4],F(0)),sum(g[2:4],F(0)),g[3],None,g[4],sum(g[4:6],F(0)),sum(g[4:7],F(0)),sum(g[4:8],F(0))]
coord_mismatch=F(0);trig_halfwidth=F(0)
for node in range(9):
 if node==4:ints=[(F(-1),F(-1)),(F(0),F(0))];qc=(F(-1),F(0))
 else:
  z=ss[node];sv=sin_iv(enc(z));cv=cos_iv(enc(z));si=(F(sv.lo,TQ),F(sv.hi,TQ));ci=(F(cv.lo,TQ),F(cv.hi,TQ));trig_halfwidth=max(trig_halfwidth,F(sv.hi-sv.lo,2*TQ),F(cv.hi-cv.lo,2*TQ))
  ints=[(-ci[1],-ci[0]),si] if node<4 else [(-ci[1],-ci[0]),(-si[1],-si[0])]
  a,b=nv[node];qc=(x[a],x[b])
 for k in range(2):coord_mismatch=max(coord_mismatch,abs(qc[k]-ints[k][0]),abs(qc[k]-ints[k][1]))
assert coord_mismatch<F(1,10**31)
# Infinity coordinate mismatch d gives chord <=2d and angle <=2*chord;
# a gap is a difference of two angles.  Thus gap error <8(d+rad).
gap_error=8*(coord_mismatch+rad)
assert gap_error<F(1,10**29)
assert trig_halfwidth<F(1,10**31)

# ---------- ninth-gap and 9D-to-8D local-ball linkage ----------
# The branch verifier uses the first eight rational center gaps above and the
# ninth center stored in the high-precision proposal file.  Certify that this
# ninth value agrees with 2*pi-sum(g[0:8]) to far better than its EPS=1e-29
# coordinate inflation.
H=json.load(open('/mnt/data/cover11_highprec_stresses.json'))
branch_center9=F(H['gaps'][8]);phi9_lo=core.S_L-sum(g,F(0));phi9_hi=core.S_U-sum(g,F(0))
branch_phi9_mismatch=max(abs(branch_center9-phi9_lo),abs(branch_center9-phi9_hi))
ninth_gap_error=branch_phi9_mismatch+8*gap_error
BRANCH_EPS=F(1,10**29)
assert branch_phi9_mismatch<F(1,10**43)
assert gap_error<BRANCH_EPS and ninth_gap_error<BRANCH_EPS
# The derivative matrix of Phi(u)=(u,2*pi-sum u) has Gram matrix I+11^T.
# Thus ||Phi(u)-Phi(u*)||_2^2-||u-u*||_2^2=(sum h_i)^2>=0.
J=[[F(i==j) for j in range(8)] for i in range(8)]+[[-F(1)]*8]
Gram=[[sum((J[k][i]*J[k][j] for k in range(9)),F(0)) for j in range(8)] for i in range(8)]
for i in range(8):
 for j in range(8):assert Gram[i][j]==F(i==j)+1

# ---------- central-to-exact derivative perturbation budget ----------
# The local central audit rounds proposal D coefficients to 110 bits.
DBITS=110;DQ=1<<DBITS;Dround_err=F(0);Dbar=[]
for idx in map(int,L['selected']):
 row=[]
 for p,q in L['D'][str(idx)]:
  z=F(int(p),int(q));v=z*DQ;n=(2*v.numerator+v.denominator)//(2*v.denominator);r=F(n,DQ);Dround_err=max(Dround_err,abs(r-z));row.append(r)
 Dbar.append(row)
assert Dround_err<F(1,10**33)
dD=Dround_err+dD_correction
assert dD<F(2,10**33)
# Pair-angle error <=8 gap_error.  Recompute the actual fixed-point
# enclosure half-width over every one of the 36 central pair sums used by the
# derivative audit (rather than borrowing a bound from the anchor angles).
central_trig_halfwidth=F(0)
for i in range(8):
 for j in range(i+1,9):
  z=sum(g[i:j],F(0));sv=sin_iv(enc(z));cv=cos_iv(enc(z))
  central_trig_halfwidth=max(central_trig_halfwidth,
      F(sv.hi-sv.lo,2*TQ),F(cv.hi-cv.lo,2*TQ))
# sin and cos are 1-Lipschitz.
trig_error=8*gap_error+central_trig_halfwidth
assert trig_error<F(1,10**28)
# The exact audit proves sum D*l < 3/2 and M3<3.07 for every row.
# Crude component and operator error bounds:
deriv_component=36*dD+F(3,2)*trig_error
grad_row_error=3*deriv_component      # sqrt(8)<3
hess_op_error=8*deriv_component
assert grad_row_error<F(1,10**26)
assert hess_op_error<F(1,10**25)

# ---------- rank/subspace and final constant transfer ----------
# Exact central audit: each row's N0 residual norm <1e-11 and sigma4(P)>=.007.
# Therefore spectral perturbation from the rank-4 projected matrix is <3e-11.
delta=F(3,10**11);sigma=F(7,1000)
assert delta<sigma/2
# Standard rank-preserving Wedin bound; use projector distance tau=1e-8.
tau=F(1,10**8)
assert 2*delta/(sigma-delta)<tau
# Unit vectors in the two subspaces may be matched within 3*tau.
vecdist=3*tau;Hnorm=F(3,2)
# sharpness: central 0.00135, gradient norm <=M3<3.07
# The central audit also proves that each central gradient's residual
# orthogonal to A0 has norm <1e-11.  Subtract it explicitly here.
central_projection_residual=F(1,10**11)
sharp=F(27,20000)-F(307,100)*vecdist-central_projection_residual-grad_row_error
assert sharp>F(13,10000)
# Hessian block transfers from central constants .0027,.0576,.5685.
quadloss=2*Hnorm*vecdist+Hnorm*vecdist*vecdist+hess_op_error
crossloss=Hnorm*(2*vecdist+vecdist*vecdist)+hess_op_error
assert F(27,10000)-quadloss>F(13,5000)
assert F(36,625)+crossloss<F(29,500)
assert F(1137,2000)+quadloss<F(57,100)
# Third derivative transfer is far below the .010096... margin.
third_extra=36*dD*F(23) # 8^(3/2)<23
assert F(3059904,1000000)+third_extra<F(307,100)

# ---------- final radius inequality ----------
kappa=F(13,10000);m=F(13,5000);cross=F(29,500);p=F(57,100);M=F(307,100);rho=F(11,5000)
end0=m/2-M*rho/6
end1=kappa-(cross+p/2)*rho-M*rho*rho/6
assert end0>0 and end1>0
# Local leaves were classified around the rational center with radial slack
# >1.4e-7, so the exact center displacement is negligible.
leaf_radial_slack=F(14,10**8) # 1.4e-7, slightly below measured exact slack
assert 3*gap_error<leaf_radial_slack

print('LOCAL TRANSFER CERTIFICATE VERIFIED')
print('Tplus minus t upper endpoint',float(TARGET-(x[TCOL]+rad)))
print('B inverse eta',float(eta),'proposal residual',float(root_res),'weight correction inf',float(dwinf))
print('gap error',float(gap_error),'ninth gap error',float(ninth_gap_error),'central trig halfwidth',float(central_trig_halfwidth),'D total error',float(dD),'grad error',float(grad_row_error),'H error',float(hess_op_error))
print('subspace tau',float(tau),'sharpness lower',float(sharp),'H transfer loss',float(quadloss))
print('local endpoint margins',float(end0),float(end1))
