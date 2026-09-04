#!/usr/bin/env python3
"""Independent set-level audit of the safe Farkas residual partition."""
from pathlib import Path
import json,pickle,sys,time
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from enumerate_orbits import canon_state,canonical_orbit
CASES=[(9,1),(9,2),(9,3),(10,1),(10,2),(11,1)]
CAND=(9,3,21015)
CAND_FACES=[(0,1,11),(0,8,10),(0,10,11),(1,2,11),(2,3,11),
            (3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),
            (6,9,10),(7,8,10),(9,10,11)]

def verify():
 t0=time.time();metric=json.load(open(ROOT/'metric_exact_cert_safe.json'))
 safe={};case_counts={}
 for B,I in CASES:
  reps=pickle.load(open(ROOT/f'triangulations_B{B}_I{I}_orbits.pkl','rb'))
  ids=list(map(int,metric['cases'][f'{B},{I}']['residual_indices']))
  assert len(ids)==len(set(ids)) and all(0<=i<len(reps) for i in ids)
  for idx in ids:safe[(B,I,idx)]=reps[idx]
  case_counts[f'{B},{I}']=len(ids)
 assert len(safe)==132 and CAND in safe
 candcanon=canonical_orbit(canon_state(CAND_FACES),9,3)
 assert candcanon==safe[CAND]
 # The exact residual certificate must contain every safe noncandidate and no other.
 exact=json.load(open(ROOT/'residual_energy_cert_safe.json'));emap={}
 for z in exact['topologies']:
  key=(int(z['B']),int(z['I']),int(z['idx']));assert key not in emap and key!=CAND
  face=canon_state(z['faces']);assert face==safe[key]
  emap[key]=face
 assert set(emap)==set(safe)-{CAND} and len(emap)==131
 report={'verified':True,'safe_residuals':len(safe),'candidate':list(CAND),
         'noncandidate_exact':len(emap),
         'case_counts':case_counts,'candidate_faces':[list(f) for f in safe[CAND]],
         'seconds':time.time()-t0}
 (ROOT/'cover12_topology_partition_verified.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':verify()
