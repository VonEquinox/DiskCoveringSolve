#!/usr/bin/env python3
"""Fail-closed, from-source replay of the n=20 computer-assisted proof.
No pre-existing success report is an input to this program.
"""
if not __debug__:
 import pathlib as _pathlib
 _base=_pathlib.Path(__file__).resolve().parent/'reports'
 for _name in ['MASTER_VERIFIED.json','MASTER_VERIFIED.json.tmp']:
  (_base/_name).unlink(missing_ok=True)
 raise RuntimeError('Assertions must be enabled: never use -O/-OO')
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import sys,os,json,subprocess,hashlib,time,shutil
ROOT=Path(__file__).resolve().parent;CORE=ROOT/'core';ENUM=ROOT/'enumeration';OUT=ROOT/'reports';RUN=ROOT/'_runtime'
OUT.mkdir(exist_ok=True);RUN.mkdir(exist_ok=True)
MASTER=OUT/'MASTER_VERIFIED.json'
CORE_SOURCES=['system20.py','verify_root20.py','interval20.py','verify_upper20.py','exact_arcs.py','anchor_isolation20.py','candidate_global20.py','compact_force20.py','domain20.py','cases20.py','verify_tree20.py','verify_interfaces20.py','verify_peeling20.py','verify_recursion_tests20.py','verify_forest20.py','test_fail_closed20.py']
INPUTS=['PROOF_zh.md','README.md','verify_all.py']+['core/'+s for s in CORE_SOURCES]+['core/'+s for s in ['root20.json','anchor_angles20.json','anchor_isolation_certificate.json','candidate20_tree.json.gz']]+['enumeration/'+s for s in ['peel_dfs20.cpp','peel_bfs20.cpp','compare_maps20.cpp','cases20_meta.json','candidate20_map.json','noncandidate20_trees.jsonl.gz']]+[f'enumeration/survivors20_B{B}.bin.gz' for B in range(12,20)]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def run(name,args):
 start=time.time();lp=OUT/(name+'.log')
 with open(lp,'w') as f:
  p=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError((name,'failed',p.returncode,'see',str(lp)))
 return {'step':name,'seconds':time.time()-start,'log_sha256':sha(lp)}
def main():
 start=time.time();MASTER.unlink(missing_ok=True);(OUT/'MASTER_VERIFIED.json.tmp').unlink(missing_ok=True)
 env=os.environ;env.pop('PYTHONOPTIMIZE',None)
 workers=int(env.get('COVER20_WORKERS','4'));assert 1<=workers<=16
 hashes={x:sha(ROOT/x) for x in INPUTS};steps=[]
 log=OUT/'FULL_REPLAY.log'
 def say(s):
  print(s,flush=True)
  with open(log,'a') as f:f.write(s+'\n')
 log.write_text('N20 from-source exact replay; no numerical optimizer is trusted.\n')
 try:
  for fn,rp in [('verify_root20.py','root20_verified.json'),('verify_upper20.py','upper20_verified.json'),('anchor_isolation20.py','anchor_isolation_verified.json'),('verify_interfaces20.py','interfaces20_verified.json')]:
   (CORE/rp).unlink(missing_ok=True);(OUT/rp).unlink(missing_ok=True)
   steps.append(run(fn,[sys.executable,'-S','-B',str(CORE/fn)]));d=read(CORE/rp);assert d['verified'] is True;shutil.copy(CORE/rp,OUT/rp);say(fn+' VERIFIED')
  cxx=env.get('CXX') or shutil.which('g++') or shutil.which('clang++');assert cxx,'C++17 compiler required'
  for name in ['peel_dfs20','peel_bfs20','compare_maps20']:
   steps.append(run('compile_'+name,[cxx,'-O3','-std=c++17',str(ENUM/(name+'.cpp')),'-o',str(RUN/name)]))
  jobs=[]
  with ThreadPoolExecutor(max_workers=workers) as ex:
   for method in ['dfs','bfs']:
    for B in range(12,20):
     path=RUN/f'{method}_B{B}.json';path.unlink(missing_ok=True)
     jobs.append(ex.submit(run,f'{method}_B{B}',[str(RUN/f'peel_{method}20'),'20',str(B),str(path)]))
   for future in as_completed(jobs):
    step=future.result();steps.append(step);say(step['step']+' exhaustive recursion completed')
  for fn,rp in [('verify_peeling20.py','peeling20_verified.json'),('verify_recursion_tests20.py','recursion_tests20_verified.json'),('verify_forest20.py','forest20_verified.json'),('test_fail_closed20.py','negative20_verified.json')]:
   (OUT/rp).unlink(missing_ok=True);steps.append(run(fn,[sys.executable,'-S','-B',str(CORE/fn)]));assert read(OUT/rp)['verified'] is True;say(fn+' VERIFIED')
  root=read(OUT/'root20_verified.json');upper=read(OUT/'upper20_verified.json');local=read(OUT/'anchor_isolation_verified.json');interfaces=read(OUT/'interfaces20_verified.json');peel=read(OUT/'peeling20_verified.json');forest=read(OUT/'forest20_verified.json');neg=read(OUT/'negative20_verified.json');reg=read(OUT/'recursion_tests20_verified.json')
  for d in [root,upper,local,interfaces,forest]:assert d['root_sha256']==hashes['core/root20.json']
  assert forest['local_certificate_sha256']==local['certificate_sha256']==hashes['core/anchor_isolation_certificate.json']
  assert root['dimension']==198 and len(root['boxes'])==2
  assert [q['rho'] for q in root['boxes']]==[str(F(1,10**60)),str(F(1,10**80))]
  assert all(0<=F(q['inclusion_ratio'])<1 and 0<=F(q['contraction'])<1 for q in root['boxes'])
  assert F(root['minimum_weight'])>0
  assert (upper['vertices'],upper['edges'],upper['faces'],upper['boundary_caps'],upper['exact_equality_center_faces'])==(33,83,51,13,19)
  assert len(upper['strict_center_faces'])==6 and F(upper['minimum_orientation'])>0 and F(upper['minimum_arc_wedge_margin'])>0
  assert F(local['anchor_Frobenius_radius'])==F(11,100) and F(local['isolation_margin'])==F(1199,25000000)>0
  assert forest['candidate_key']==[interfaces['candidate_B'],interfaces['candidate_index']]
  total=peel['surviving_orbits'];assert total==forest['surviving_orbits'] and total>0
  nc=forest['noncandidate_forest'];cand=forest['candidate_tree'];assert nc['cases']+1==total
  assert nc['certificate_sha256']==hashes['enumeration/noncandidate20_trees.jsonl.gz'] and cand['certificate_sha256']==hashes['core/candidate20_tree.json.gz']
  assert F(nc['minimum_force_margin'])>0 and F(cand['minimum_force_margin'])>0 and F(cand['minimum_local_margin'])>0
  assert cand['counts']['LOCAL']>0 and nc['counts']['LOCAL']==cand['counts']['EMPTY']==nc['counts']['EMPTY']==0
  assert peel['unresolved_combinatorial_branches']==forest['unresolved_leaves']==0
  assert reg['pruning_disabled'] is True and reg['cases']>=29 and neg['negative_test_count']>=30
  assert hashes=={x:sha(ROOT/x) for x in INPUTS},'An accepting input changed during replay'
  out={'verified':True,'theorem':'r_20 = sqrt(t_star), R_20 = 1/sqrt(t_star)',
       'global_optimality_certificate_chain_closed':True,'kkt_dimension':198,'decimal_intervals':root['decimal_intervals'],
       'root_sha256':hashes['core/root20.json'],'covering_upper_verified':True,
       'upper_cover':{k:upper[k] for k in ['vertices','edges','faces','boundary_caps','exact_equality_center_faces']},
       'enumeration_method':'Exhaustive root-face peeling with safe existing-edge angular cuts; independent root/region-order/BFS/Floyd replay; full canonical survivor-set comparison.',
       'exhaustive_partial_topology_search_closed':True,'surviving_orbits':total,'rooted_survivors':peel['rooted_survivors'],'families':peel['families'],
       'noncandidate_cases':nc['cases'],'candidate_cases':1,'noncandidate_forest':nc,'candidate_tree':cand,
       'unresolved_combinatorial_branches':0,'unresolved_leaves':0,'anchor_isolation_radius':local['anchor_Frobenius_radius'],
       'negative_tests':neg,'recurrence_regression_cases':reg['cases'],'input_sha256':hashes,'steps':steps,
       'seconds':time.time()-start,'completed_utc':datetime.now(timezone.utc).isoformat(),
       'trust_boundary':'Written convex geometry, planar-cell topology and exhaustive peeling induction; Python integer/Fraction arithmetic, C++17 compiler. Discovery optimizers are not trusted. Not a Lean/Coq formalization or third-party peer review.'}
  tmp=Path(str(MASTER)+'.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(MASTER)
  say('R20 COMPLETE CERTIFICATE CHAIN VERIFIED')
  say(f"{nc['cases']} + 1 = {total} surviving orbits after exhaustive partial cuts; unresolved leaves = 0")
  return out
 except BaseException:
  MASTER.unlink(missing_ok=True);Path(str(MASTER)+'.tmp').unlink(missing_ok=True)
  raise
if __name__=='__main__':main()
