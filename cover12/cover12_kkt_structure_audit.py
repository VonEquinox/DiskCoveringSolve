#!/usr/bin/env python3
"""Symbolic audit connecting the 14-variable candidate system to the full
D3-symmetric contact graph, free-node equilibrium, and anchor stationarity.

All identities are polynomial identities over Q.  No numerical tolerance is
used.  The Groebner bases below only certify ideal membership of the displayed
contact/equilibrium polynomials.
"""
from __future__ import annotations
from pathlib import Path
import json,time,sys
import sympy as sp
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT));T0=time.time()
import cover12_symmetric_kkt_exact as KKT

r,u,v,w,c,s,z,P,Q,S,A,B,C,D=sp.symbols('r u v w c s z P Q S A B C D')
VARS=[r,u,v,w,c,s,z,P,Q,S,A,B,C,D]
h=sp.Rational(1,2)
L=(c*c-s*s)*h+z*c*s
M=z*(c*c-s*s)*h-c*s
E=[
 c*c+s*s-1,
 u*u-2*u*c+1-r*r,
 v*v-2*v*L+1-r*r,
 w-u*u+r*r+r,
 v*v-(2*r+w)*v+w*(w+r),
 z*z-3,
 2*A*(u*c-(w+r))-B*r,
 C*r-2*D*(v-r-w*h),
 2*P*(L-v)-C*r,
 Q-S-A*(w+r),
 (S+Q)*(c-u)+A*((w+r)*c-u),
 B*r-D*(2*w-v+r),
 Q*u*s-P*v*M,
 6*(P+Q+S+A+D)+3*(B+C)-1,
]

def sympy_from_sparse(poly):
 out=0
 for aa,ex in poly:
  q=sp.Rational(aa.numerator,aa.denominator)
  for x,k in zip(VARS,ex):q*=x**k
  out+=q
 return sp.expand(out)
assert len(KKT.POLYS)==len(E)==14
assert all(sp.expand(sympy_from_sparse(a)-b)==0 for a,b in zip(KKT.POLYS,E))

def add(a,b):return (sp.expand(a[0]+b[0]),sp.expand(a[1]+b[1]))
def sub(a,b):return (sp.expand(a[0]-b[0]),sp.expand(a[1]-b[1]))
def mul(q,a):return (sp.expand(q*a[0]),sp.expand(q*a[1]))
def dot(a,b):return sp.expand(a[0]*b[0]+a[1]*b[1])
def rot120(a):return (sp.expand(-(a[0]+z*a[1])*h),sp.expand((z*a[0]-a[1])*h))
def rot240(a):return (sp.expand((-a[0]+z*a[1])*h),sp.expand(-(z*a[0]+a[1])*h))
def reflect60(a):return (sp.expand((-a[0]+z*a[1])*h),sp.expand((z*a[0]+a[1])*h))

gp=(u*c,u*s);gm=(u*c,-u*s)
rawC=[gp,gm,rot120(gp),rot120(gm),rot240(gp),rot240(gm),
      (-v,0),rot120((-v,0)),rot240((-v,0)),
      (w,0),rot120((w,0)),rot240((w,0))]
CENTERS=[rawC[i] for i in [6,5,4,7,1,0,8,3,2,9,10,11]]
q0=(sp.Integer(1),sp.Integer(0));q1=(c*c-s*s,2*c*s);q2=reflect60(q1)
rawQ=[q0,q1,q2,rot120(q0),rot120(q1),rot120(q2),
      rot240(q0),rot240(q1),rot240(q2)]
ANCHORS=[rawQ[i] for i in [5,6,7,8,0,1,2,3,4]]
FACES=[(0,1,11),(0,8,10),(0,10,11),(1,2,11),(2,3,11),
       (3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),
       (6,9,10),(7,8,10),(9,10,11)]
ACTIVE=[(1,2,11),(4,5,9),(7,8,10),
        (0,10,11),(3,9,11),(6,9,10)]
pgen=(w+r,0);pax=(-v+r,0)
WITNESS=[rot240(pgen),pgen,rot120(pgen),pax,rot120(pax),rot240(pax)]

# Geometry ideal: the first six equations.  It proves all unit-circle and
# active-contact equations for all D3 copies, not merely representatives.
Ggeom=sp.groebner(E[:6],r,u,v,w,c,s,z,order='grevlex',domain=sp.QQ)
def rem0(G,x):
 rr=G.reduce(sp.expand(x))[1]
 assert rr==0,sp.factor(rr)
trig_identity=sp.expand(L*L+M*M-1);rem0(Ggeom,trig_identity)
unit=[]
for q in ANCHORS:
 x=dot(q,q)-1;rem0(Ggeom,x);unit.append(x)
boundary=[]
for i in range(9):
 for j in (i,(i+1)%9):
  x=dot(sub(ANCHORS[i],CENTERS[j]),sub(ANCHORS[i],CENTERS[j]))-r*r
  rem0(Ggeom,x);boundary.append(x)
face_contacts=[]
for k,fa in enumerate(ACTIVE):
 for j in fa:
  x=dot(sub(WITNESS[k],CENTERS[j]),sub(WITNESS[k],CENTERS[j]))-r*r
  rem0(Ggeom,x);face_contacts.append(x)

# Build the exact 36-edge active graph and its force-balance polynomials.
weights={'P':P,'Q':Q,'S':S,'A':A,'B':B,'C':C,'D':D}
pat=[('P','Q'),('S','S'),('Q','P')]*3
fc=[(0,0) for _ in range(12)];fw=[(0,0) for _ in range(6)]
for i,(left,right) in enumerate(pat):
 for j,name in ((i,left),((i+1)%9,right)):
  fc[j]=add(fc[j],mul(weights[name],sub(ANCHORS[i],CENTERS[j])))
for k,fa in enumerate(ACTIVE):
 for j in fa:
  name=('B' if j>=9 else 'A') if k<3 else ('D' if j>=9 else 'C')
  fc[j]=add(fc[j],mul(weights[name],sub(WITNESS[k],CENTERS[j])))
  fw[k]=add(fw[k],mul(weights[name],sub(CENTERS[j],WITNESS[k])))

# Each D3 orbit is reduced against precisely the equations assigned to that
# node type.  This avoids trusting a hand-written force interpretation.
def audit_group(polys,ideal):
 G=sp.groebner(ideal,*VARS,order='grevlex',domain=sp.QQ)
 for x in polys:rem0(G,x)
 return len(polys)
counts={}
counts['generic_center_components']=audit_group([x for i in [1,2,4,5,7,8] for x in fc[i]],[E[0],E[5],E[9],E[10]])
counts['axial_center_components']=audit_group([x for i in [0,3,6] for x in fc[i]],[E[0],E[5],E[8]])
counts['inner_center_components']=audit_group([x for i in [9,10,11] for x in fc[i]],[E[5],E[11]])
counts['generic_witness_components']=audit_group([x for i in [0,1,2] for x in fw[i]],[E[5],E[6]])
counts['axial_witness_components']=audit_group([x for i in [3,4,5] for x in fw[i]],[E[5],E[7]])

# Tangential force at every unit-circle anchor.  Vanishing is exactly the
# angular gradient condition for the Dirichlet energy after free-node balance.
tan=[]
for i,(left,right) in enumerate(pat):
 f=add(mul(weights[left],sub(ANCHORS[i],CENTERS[i])),
       mul(weights[right],sub(ANCHORS[i],CENTERS[(i+1)%9])))
 q=ANCHORS[i];tan.append(dot(f,(-q[1],q[0])))
counts['anchor_tangential_components']=audit_group(tan,[E[0],E[5],E[12]])

# Exact total edge weight; together with all 36 contact identities it gives
# sum_e w_e |x_e-y_e|^2 = r^2.
total=6*(P+Q+S)+3*(2*A+B)+3*(C+2*D)
assert sp.expand(total-1-E[13])==0
assert len(boundary)+len(face_contacts)==36

report={
 'verified':True,
 'kkt_polynomials_match':True,'secondary_trig_identity':True,
 'unit_anchor_identities':len(unit),
 'active_boundary_contacts':len(boundary),
 'active_face_contacts':len(face_contacts),
 'active_edges':36,
 'free_nodes':18,
 'force_components':sum(counts.values())-counts['anchor_tangential_components'],
 'anchor_tangential_components':counts['anchor_tangential_components'],
 'normalization_identity':True,
 'energy_identity':'sum_e w_e |x_e-y_e|^2 = r^2',
 'faces':[list(f) for f in FACES],
 'active_faces':[list(f) for f in ACTIVE],
 'seconds':time.time()-T0,
 'component_counts':counts,
}
(ROOT/'cover12_kkt_structure_audit_verified.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
