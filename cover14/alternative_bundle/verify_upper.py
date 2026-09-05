#!/usr/bin/env python3
"""Exact continuous coverage of the unit disk by the algebraically isolated n=14 candidate.
Six extra tight faces follow from verified symmetries and elementary exact squared-distance identities.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter,defaultdict,deque
import json,time,hashlib
from framework import *
from intervals import *
from verify_symmetry import boxes,CX,CT,PX,PT
if not __debug__:raise RuntimeError('Assertions must be enabled')
BASE=Path(__file__).resolve().parent

def verify():
 st=time.time();cert,x,rho,X,pos=boxes();Q=pos[:B];C=pos[B:B+N_DISKS];P=pos[B+N_DISKS:];t=X[TID]
 # verify_symmetry proves this exact value, not merely an interval containing it.
 Q[5]=pt(-1,0);edge_set={tuple(sorted(e)) for e in EDGES}
 for i in range(B):
  for c in (i,(i+1)%B):assert tuple(sorted((i,B+c))) in edge_set
 for k,fa in enumerate(ACTIVE_FACES):
  for c in fa:assert tuple(sorted((B+N_DISKS+k,B+c))) in edge_set
 order=[];cone=[]
 quad=[(1,0),(1,1),(1,1),(-1,1),(-1,1),(-1,0),(-1,-1),(-1,-1),(1,-1),(1,-1)]
 for p,(sx,sy) in zip(Q,quad):
  assert (sx*p[0]).lo>0
  if sy:assert (sy*p[1]).lo>0
  else:assert p[1]==I.point(0)
 for i in range(B):
  cp=cross(Q[i],Q[(i+1)%B]);dp=dot(Q[i],Q[(i+1)%B]);assert cp.lo>0 and dp.lo>0;order.extend([cp.lo,dp.lo])
  a=cross(Q[(i-1)%B],C[i]);b=cross(C[i],Q[i]);assert a.lo>0 and b.lo>0;cone.extend([a.lo,b.lo])
 # Positive-stress graph incidences needed for the six tight-face identities.
 def neighbors(v):return sorted((w if u==v else u) for u,w in EDGES if v in (u,w))
 assert neighbors(B+3)==[2,3]
 assert neighbors(B+11)==[B+N_DISKS+2,B+N_DISKS+4]
 # Vertical reflection is half_turn o reflection_x.
 CV=[CT[CX[j]] for j in range(14)]
 PV=[PT[PX[j]] for j in range(NN)]
 assert PV[2]==3 and PV[B+3]==B+3 and PV[B+11]==B+11
 assert PV[B+N_DISKS+2]==B+N_DISKS+4
 assert PX[B+10]==B+10 and PX[B+11]==B+13
 assert Q[2][0].lo>0 and P[2][0].lo>0 and (C[3][1]-C[11][1]).lo>0
 assert ACTIVE_FACES[2]==(2,10,11) and ACTIVE_FACES[4]==(4,11,12)
 # Exact implications (also proved in PROOF_zh.md):
 # C3=(0,z), Q2=(r,z), C11=(0,b), P2=(r,b), C10=(a,0), C13=(0,-b).
 # The two y equalities follow from the above two-node equilibrium equations and positive weights.
 # Both x coordinates Q2_x and P2_x are positive with square t, hence are the same r.
 # For C2=(x,y), V=(x-r,y) is a common point of disks 2,3,11:
 # |V-C2|^2=t, |V-C3|^2=|C2-Q2|^2=t, |V-C11|^2=|C2-P2|^2=t.
 # W=(a-r,0) is common to 10,11,13:
 # |W-C10|^2=t, |W-C11|^2=|W-C13|^2=|C10-P2|^2=t.
 # Each distance on the right is a declared equation in the root framework.
 for e in [(2,B+2),(B+N_DISKS+2,B+2),(B+N_DISKS+2,B+10),(2,B+3),(B+N_DISKS+2,B+11)]:assert tuple(sorted(e)) in edge_set
 group=[list(range(14)),CX,CT,CV]
 extra={tuple(sorted(g[j] for j in f)) for g in group for f in [(2,3,11),(10,11,13)]}
 assert extra=={ALL_FACES[k] for k in [3,5,10,12,14,15]}
 slack={(0,1,10):pt(F(685075,1000000),0),(5,6,12):pt(-F(685075,1000000),0)}
 assert set(ALL_FACES)==set(ACTIVE_FACES)|extra|set(slack)
 assert len(set(ACTIVE_FACES))+len(extra)+len(slack)==len(ALL_FACES)==16
 sm=[]
 for fa,p in slack.items():
  for c in fa:
   m=t.lo-dist2(p,C[c]).hi;assert m>0;sm.append(m)
 # An oriented triangulated disk filling the anchor decagon.
 faces=[]
 for i in range(B):
  # For i=3,8 the radial triangle is exactly degenerate: center is the chord midpoint.
  if i not in (3,8):faces.append(((i-1)%B,i,B+i))
  faces.append((i,B+(i+1)%B,B+i))
 for fa in ALL_FACES:
  ids=[B+j for j in fa];v=orient(C[fa[0]],C[fa[1]],C[fa[2]])
  if v.hi<0:ids[1],ids[2]=ids[2],ids[1]
  else:assert v.lo>0
  faces.append(tuple(ids))
 def V(i):return Q[i] if i<B else C[i-B]
 om=[];ec=Counter();directed=Counter()
 for fa in faces:
  o=orient(*(V(v) for v in fa));assert o.lo>0;om.append(o.lo)
  for a,b in zip(fa,fa[1:]+fa[:1]):ec[tuple(sorted((a,b)))]+=1;directed[a,b]+=1
 boundary=[0,1,2,B+3,3,4,5,6,7,B+8,8,9]
 outer_oriented=list(zip(boundary,boundary[1:]+boundary[:1]));outer={tuple(sorted(e)) for e in outer_oriented}
 assert {e for e,n in ec.items() if n==1}==outer and all(n in (1,2) for n in ec.values())
 for (a,b),n in ec.items():
  if n==2:assert directed[a,b]==directed[b,a]==1
 for i,j in outer_oriented:assert directed[i,j]==1 and directed[j,i]==0
 assert len(faces)==34 and len(ec)==57 and 24-len(ec)+len(faces)==1
 # Links certify a connected combinatorial disk with the asserted boundary cycle.
 adj=defaultdict(set)
 for a,b in ec:adj[a].add(b);adj[b].add(a)
 def connected(g):
  s={next(iter(g))};todo=list(s)
  for u in todo:
   for v in g[u]:
    if v not in s:s.add(v);todo.append(v)
  return s==set(g)
 assert connected(adj) and set(adj)==set(range(24))
 for v in range(24):
  g=defaultdict(set)
  for fa in faces:
   if v in fa:
    a,b=[u for u in fa if u!=v];g[a].add(b);g[b].add(a)
  assert connected(g);ds=[len(a) for a in g.values()]
  if v in boundary:assert ds.count(1)==2 and all(d in (1,2) for d in ds)
  else:assert all(d==2 for d in ds)
 out={'verified':True,'vertices':24,'edges':57,'triangles':34,'boundary_vertices':12,'omitted_degenerate_radial_faces':2,'center_faces':16,'positive_stress_faces':8,'extra_tight_faces':6,'strict_slack_faces':2,'boundary_order_margin':str(min(order)),'center_cone_margin':str(min(cone)),'triangle_orientation_margin':str(min(om)),'slack_witness_margin':str(min(sm)),'root_sha256':hashlib.sha256((BASE/'root_certificate.json').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (BASE/'upper_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':verify()
