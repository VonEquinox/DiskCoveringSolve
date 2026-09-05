#!/usr/bin/env python3
"""Fresh end-to-end acceptance. No discovery solver or stored success flag is trusted."""
from pathlib import Path
import sys,json,subprocess,time,hashlib,datetime
R=Path(__file__).resolve().parent;REPORT=R/'reports'/'MASTER_VERIFIED.json'
CORE=['geometry19.py','exact_arcs.py','exact_force19.py','verify_interfaces19.py','verify_candidate19.py','verify_enumeration19.py','verify_forest19.py','test_fail_closed19.py']
INPUTS=['PROOF_zh.md','README.md','verify_all.py']+['core/'+p for p in CORE]+['enumeration/peel_enum.cpp','enumeration/peel_audit.cpp','enumeration/metric19_residuals.json','enumeration/candidate19_map.json','enumeration/noncandidate19_trees.jsonl.gz']+[f'enumeration/B{b}_I{19-b}.txt' for b in range(12,19)]
def hashes():return {p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in INPUTS}
def main():
 REPORT.parent.mkdir(exist_ok=True);REPORT.unlink(missing_ok=True)
 if not __debug__:raise RuntimeError('Do not disable assertions')
 if len(sys.argv)!=1:raise RuntimeError('This verifier has no skip or partial mode')
 start=time.time();before=hashes();steps=[]
 def step(program,output):
  t=time.time();p=subprocess.run([sys.executable,'-S','-B',str(R/'core'/program)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,cwd=R)
  log=R/'reports'/f'{program}.log';log.write_text(p.stdout,encoding='utf8')
  if p.returncode:raise RuntimeError(f'{program} failed:\n{p.stdout}')
  d=json.loads((R/'core'/output).read_text());assert d.get('verified') is True
  steps.append({'step':program,'seconds':time.time()-t,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()});print('PASSED',program,flush=True);return d
 try:
  interface=step('verify_interfaces19.py','interfaces19_verified.json')
  candidate=step('verify_candidate19.py','candidate19_verified.json')
  enumeration=step('verify_enumeration19.py','enumeration19_verified.json')
  forest=step('verify_forest19.py','forest19_verified.json')
  negative=step('test_fail_closed19.py','negative19_verified.json')
  assert interface['radius_squared']==candidate['exact_radius_squared']=='1/13'
  assert candidate['global_strict_convexity_verified'] and candidate['positive_rods']==42 and candidate['boundary_anchors']==12
  assert enumeration['N']==19 and [v['B'] for v in enumeration['families']]==list(range(12,19))
  assert enumeration['surviving_orbits']==forest['noncandidate_cases']+forest['candidate_cases'] and forest['candidate_cases']==1
  assert forest['residual_sha256']==enumeration['residual_sha256']==before['enumeration/metric19_residuals.json']
  assert forest['unresolved_leaves']==0 and candidate['geometry_sha256']==before['core/geometry19.py']
  assert hashes()==before,'An accepting input changed during the replay'
  out={'verified':True,'theorem':'r_19 = 1/sqrt(13), R_19 = sqrt(13)','global_optimality_certificate_chain_closed':True,'exact_squared_radius':'1/13','r_strict_interval':candidate['r_interval'],'R_strict_interval':candidate['R_interval'],'covering_upper_verified':True,'upper_cover':{k:candidate[k] for k in ['full_centers','boundary_anchors','center_faces','covering_chain_faces','positive_triangles','degenerate_triangles','boundary_caps']},'candidate_global_lower':{k:candidate[k] for k in ['positive_graph_vertices','positive_rods','radial_multiplier','minimum_ldl_pivot','circulant_lower_weights','nonconstant_spectral_gap_lower']},'surviving_topology_orbits':enumeration['surviving_orbits'],'primary_partial_states':enumeration['primary_states'],'independent_partial_states':enumeration['independent_states'],'enumeration_families':enumeration['families'],'unpruned_recurrence_tests':len(enumeration['unpruned_count_regressions']),'noncandidate_forest':forest,'candidate_cases':1,'unresolved_combinatorial_branches':0,'unresolved_leaves':0,'rejection_tests':negative,'accepting_inputs_sha256':before,'steps':steps,'seconds':time.time()-start,'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'trust_boundary':'Written convex-geometric, topological and finite-dimensional analytic reductions in PROOF_zh.md; Python standard-library integer/Fraction arithmetic; a C++17 compiler. No floating-point optimizer, precomputed success flag, or previous r_n optimality theorem is trusted.','external_peer_review_completed':False,'proof_assistant_formalization_completed':False}
  REPORT.write_text(json.dumps(out,indent=2),encoding='utf8');print('R19 COMPLETE CERTIFICATE CHAIN VERIFIED',flush=True);print(json.dumps({k:out[k] for k in ['theorem','r_strict_interval','surviving_topology_orbits','primary_partial_states','independent_partial_states','unresolved_leaves','seconds']},indent=2),flush=True)
  return out
 except BaseException:
  REPORT.unlink(missing_ok=True);raise
if __name__=='__main__':main()
