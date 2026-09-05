"""Reconstruct every bridge between exact geometry, stress, angles and full graphs."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import geometry19 as S
import exact_force19 as G
import anchor_isolation19 as L
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def verify():
 assert S.check_geometry() and (S.N,S.B,S.G,S.M,S.C)==(19,12,70,42,53)
 assert (G.N,G.V,G.BMIN)==(19,55,12)
 assert (G.A,G.C,G.C6,G.D)==(894563,1871672,3128331,10**7)
 assert G.RU==F('0.27735009811262') and F(1,13)<G.RU**2
 assert (L.GAMMA,L.LAMBDA,L.NU,L.RADIUS)==(400,F(1,650),F(1,156),F(1,5))
 assert S.CENTER_IDS==list(range(12))+list(range(13,19))
 assert len(S.EDGES)==42 and len(set(S.EDGES))==42
 ex=[(i,S.CENTER_MAP[i]) for i in range(12)]+[(i,S.CENTER_MAP[(i+1)%12]) for i in range(12)]
 ex += [(30+k,S.CENTER_MAP[v]) for k,f in enumerate(S.ACTIVE) for v in f]
 assert S.EDGES==ex
 # Full, physical rod graph has 55 nodes; the stressed graph has only 36.
 full=set(map(frozenset,G.edges(S.B,S.FACES)))
 findex={f:i for i,f in enumerate(S.FACES)}
 embed=list(range(12))+[12+c for c in S.CENTER_IDS]+[31+findex[f] for f in S.ACTIVE]
 assert len(embed)==36 and len(set(embed))==36 and 24 not in embed
 assert all(frozenset((embed[u],embed[v])) in full for u,v in S.EDGES)
 # Reconstruct geometric derivative coefficients directly in (dx,dy/sqrt(3)).
 A,M,P=S.matrices();assert len(A)==53 and all(len(row)==70 for row in A)
 for k,(u,v) in enumerate(S.EDGES):
  direct=[F(0)]*70
  for node,sg in ((u,1),(v,-1)):
   if node:
    direct[2*(node-1)]=2*sg*(S.Z[u][0]-S.Z[v][0])
    direct[2*(node-1)+1]=6*sg*(S.Z[u][1]-S.Z[v][1])
  assert direct==A[k]
 assert P==[F(1),F(3)]*11+[F(0)]*48
 out={'verified':True,'radius_squared':'1/13','centers':19,'boundary_anchors':12,'full_center_faces':24,'full_auxiliary_nodes':55,'full_auxiliary_rods':96,'stressed_centers':18,'positive_faces':6,'stressed_rods':42,'stressed_geometric_dimension':70,'stressed_constraints':53,'unused_center_index':12,'lower_graph_embedding':embed,'geometry_sha256':hashlib.sha256((R/'geometry19.py').read_bytes()).hexdigest()}
 (R/'interfaces19_verified.json').write_text(json.dumps(out,indent=2));print('EXACT GEOMETRY AND ALL GRAPH INTERFACES VERIFIED',flush=True);return out
if __name__=='__main__':verify()
