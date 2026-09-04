#!/usr/bin/env python3
"""Exact cross-checks linking residual orbit indices to their special proofs."""
import json,pickle

ROOT='/mnt/data'
def faceset(fs):return {tuple(sorted(map(int,f))) for f in fs}

metric=json.load(open(f'{ROOT}/cover11_metric_Tplus.json'))
reps92=pickle.load(open(f'{ROOT}/triangulations_B9_I2.pkl','rb'))
reps101=pickle.load(open(f'{ROOT}/triangulations_B10_I1.pkl','rb'))
residual={}
for key,block in metric['cases'].items():
 B,I=map(int,key.split(','))
 for z in block['residual']:
  idx=int(z['idx']);fs=faceset(z['faces'])
  residual[(B,I,idx)]=fs
  reps=reps92 if (B,I)==(9,2) else reps101 if (B,I)==(10,1) else None
  if reps is not None:assert fs==faceset(reps[idx]),(B,I,idx,'metric residual differs from stored orbit')

# Ten-boundary special orbit: exactly the wheel with center vertex 10.
wheel_key=(10,1,1003);wheel={tuple(sorted((i,(i+1)%10,10))) for i in range(10)}
assert wheel_key in residual and residual[wheel_key]==wheel
assert faceset(reps101[1003])==wheel

# Candidate full core triangulation.  Boundary order is the certified center
# cycle [1,3,5,7,9,10,8,6,4]; centers 0 and 2 are the two interior vertices.
boundary_centers=[1,3,5,7,9,10,8,6,4]
boundary_index={c:i for i,c in enumerate(boundary_centers)}
active=[(0,1,3),(0,1,4),(0,2,5),(0,2,6),(0,3,5),(0,4,6),(2,5,7),(2,6,8),(2,9,10)]
slack=[(2,7,9),(2,8,10)]
def abstract_vertex(c):
 if c==0:return 9
 if c==2:return 10
 return boundary_index[c]
candidate=faceset(tuple(abstract_vertex(c) for c in f) for f in active+slack)
candidate_key=(9,2,2547);assert candidate_key in residual
orbit=faceset(reps92[2547]);assert orbit==residual[candidate_key]
# Explicit witness requested by the paper audit: map residual boundary vertex i
# to i+3 mod 9 and swap the two interior vertices 9<->10.
def witness(v):return (v+3)%9 if v<9 else 19-v
assert faceset(tuple(witness(v) for v in f) for f in orbit)==candidate

# Removing precisely the candidate and wheel leaves the 53 ordinary residuals.
assert len(residual)==55 and len(set(residual)-{candidate_key,wheel_key})==53
report={
 'status':'passed',
 'candidate':{'B':9,'I':2,'orbit':2547,'boundary_shift':3,'boundary_reflection':False,'interior_swap':True},
 'wheel':{'B':10,'I':1,'orbit':1003,'exact_wheel':True},
 'residual_total':55,'ordinary_residuals':53,
}
open(f'{ROOT}/cover11_topology_linkage_verify.json','w').write(json.dumps(report,indent=2))
print('SPECIAL TOPOLOGY LINKAGE VERIFIED')
print('candidate orbit (9,2,2547): shift=3, interior swap=True')
print('wheel orbit (10,1,1003): exact 10-wheel')
