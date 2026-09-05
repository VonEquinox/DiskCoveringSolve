#!/usr/bin/env python3
"""Exact root-box transport proving two candidate symmetries; no imposed symmetry in the KKT solve.
The omitted q0 multiplier is recovered using the zero total torque identity.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time,hashlib
from framework import *
from verify_root import load
from intervals import *
if not __debug__:raise RuntimeError('Assertions must be enabled')
BASE=Path(__file__).resolve().parent

def permutation(cperm,qperm):
 assert sorted(cperm)==list(range(14)) and sorted(qperm)==list(range(10))
 fa={tuple(sorted(f)):j for j,f in enumerate(ACTIVE_FACES)}
 pperm=[fa[tuple(sorted(cperm[c] for c in f))] for f in ACTIVE_FACES]
 perm=qperm+[B+c for c in cperm]+[B+N_DISKS+p for p in pperm]
 assert sorted(perm)==list(range(NN))
 lookup={tuple(sorted(e)):k for k,e in enumerate(EDGES)}
 ep=[lookup[tuple(sorted((perm[u],perm[v])))] for u,v in EDGES]
 assert sorted(ep)==list(range(NE))
 return perm,ep
CX=[(1-j)%10 for j in range(10)]+[10,13,12,11]
QX=[(-j)%10 for j in range(10)]
CT=[(j+5)%10 for j in range(10)]+[12,13,10,11]
QT=[(j+5)%10 for j in range(10)]
PX,EX=permutation(CX,QX);PT,ET=permutation(CT,QT)

def boxes():
 cert,x,Y,rho=load();X=[I(v-rho,v+rho) for v in x];pos=[pt(1,0)]+[(X[2*i],X[2*i+1]) for i in range(NN-1)]
 return cert,x,rho,X,pos

def verify():
 st=time.time();cert,x,rho,X,pos=boxes();big=F(1,10**50)
 assert len(X)==116 and rho<big
 # Half-gradient resultant at fixed anchor; all other torques vanish at a KKT root.
 # sum_nodes cross(z, sum_edges w(z-z_other))=0 identically, so this resultant is radial.
 mu0=-sum((X[WSTART+k]*(1-pos[(v if u==0 else u)][0]) for k,(u,v) in enumerate(EDGES) if u==0 or v==0),I.point(0))
 mus=[mu0]+X[MUSTART:]
 ratios={}
 def transport(perm,ep,sign,name):
  pp=[tuple(sign[d]*pos[perm[i]][d] for d in range(2)) for i in range(NN)]
  assert pp[0]==pt(1,0)
  xx=[v for p in pp[1:] for v in p]+[X[TID]]+[X[WSTART+ep[k]] for k in range(NE)]+[mus[perm[i]] for i in range(1,B)]
  assert len(xx)==N
  err=max(max(abs(v.lo-x[i]),abs(v.hi-x[i])) for i,v in enumerate(xx))
  assert err<big
  ratios[name]=str(err/big)
 # Reflection fixes q0; the interval image is in the larger uniqueness box.
 transport(PX,EX,(1,-1),'reflection_x')
 assert PX[5]==5 and pos[5][0].hi<0
 # Reflection gives q5_y=0; its unit-circle equation and sign give q5=(-1,0).
 pos[5]=pt(-1,0)
 # Now the half-turn + relabeling also fixes q0 exactly.
 transport(PT,ET,(-1,-1),'half_turn')
 out={'verified':True,'reflection_x_node_permutation':PX,'half_turn_node_permutation':PT,'reflection_x_rod_permutation':EX,'half_turn_rod_permutation':ET,'large_box_radius':str(big),'transport_relative_error_bounds':ratios,'root_sha256':hashlib.sha256((BASE/'root_certificate.json').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (BASE/'symmetry_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return out
if __name__=='__main__':verify()
