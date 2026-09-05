#!/usr/bin/env python3
"""Fresh fail-closed replay, including regenerated enumeration and independent BFS audit."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CORE=ROOT/'core';ENUM=ROOT/'enumeration';REPORT=ROOT/'reports';RUNTIME=ROOT/'_runtime';MASTER=REPORT/'MASTER_VERIFIED.json'
if not __debug__:
 MASTER.unlink(missing_ok=True)
 raise RuntimeError('Run without -O/-OO/PYTHONOPTIMIZE')
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from math import comb
import hashlib,json,os,shutil,subprocess,sys,time
INPUTS=['PROOF_zh.md','README.md','verify_all.py',
 *['core/'+n for n in ['system15.py','root15.json','verify_root15.py','interval15.py','exact_arcs.py','verify_upper15.py','anchor_angles15.json','anchor_isolation15.py','anchor_isolation_certificate.json','exact_force15.py','candidate_global15.py','candidate15_tree.json.gz','verify_forest15.py','verify_partition15.py','test_fail_closed15.py']],
 *['enumeration/'+n for n in ['enumerate15.cpp','state_store.hpp','audit15.cpp','metric15_certificate.json','metric15_residuals.json','candidate15_map.json','noncandidate15_trees.jsonl.gz']]]

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
@lru_cache(None)
def labelled(B,I):
 if B==2:return int(I==0)
 ans=sum(comb(I,j)*labelled(k,j)*labelled(B-k+1,I-j) for k in range(2,B) for j in range(I+1))
 if I:ans+=I*(labelled(B+1,I-1)-sum(comb(I-1,j)*labelled(3,j)*labelled(B,I-1-j) for j in range(I)))
 return ans

def run_python(name,outputs):
 for n in outputs:
  (CORE/n).unlink(missing_ok=True);(REPORT/n).unlink(missing_ok=True)
 log=REPORT/(name+'.log');st=time.time()
 with log.open('wb') as f:subprocess.run([sys.executable,'-S','-B',str(CORE/name)],cwd=CORE,stdout=f,stderr=subprocess.STDOUT,check=True)
 result={}
 for n in outputs:
  p=CORE/n;d=json.loads(p.read_text());assert d['verified'] is True;result[n]=d;shutil.move(p,REPORT/n)
 print('PASSED',name,flush=True)
 return result,{'step':name,'seconds':time.time()-st,'log_sha256':sha(log)}

def family(B):
 st=time.time();fam=f'B{B}_I{15-B}';prefix=RUNTIME/fam;out=REPORT/f'audit15_B{B}';log=REPORT/(fam+'.log')
 for p in [prefix.with_suffix('.masks'),prefix.with_suffix('.json'),out.with_suffix('.json'),out.with_suffix('.residuals')]:p.unlink(missing_ok=True)
 c=json.load(open(ENUM/'metric15_certificate.json'))[fam];cs=c['certified_signatures'];lines=[f'10306 21954 100000 1000000000000 {len(cs)}']
 for key,nums in cs.items():
  masks=list(map(int,key.split(','))) if key else [];lines.extend([str(len(masks)),' '.join(map(str,masks)),' '.join(map(str,nums))])
 cuts=prefix.with_suffix('.farkas');cuts.write_text('\n'.join(lines)+'\n')
 with log.open('wb') as f:
  subprocess.run([str(RUNTIME/'enumerate15'),str(B),str(15-B),str(prefix)],stdout=f,stderr=subprocess.STDOUT,check=True)
  assert sha(prefix.with_suffix('.masks'))==c['source_masks_sha256'],'regenerated enumeration differs'
  subprocess.run([str(RUNTIME/'audit15'),str(prefix.with_suffix('.masks')),str(cuts),str(labelled(B,15-B)),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True)
 d=json.load(open(out.with_suffix('.json')));assert d['verified'] is True and int(d['orbit_mass'])==labelled(B,15-B)
 print('PASSED regenerated enumeration and independent BFS audit',B,15-B,'orbits',d['representatives'],flush=True)
 return {'step':fam,'seconds':time.time()-st,'log_sha256':sha(log),'labelled_count':labelled(B,15-B)}

def main():
 REPORT.mkdir(exist_ok=True);RUNTIME.mkdir(exist_ok=True);MASTER.unlink(missing_ok=True);st=time.time();steps=[]
 if sys.byteorder!='little':raise RuntimeError('The supplied binary interchange format requires a little-endian machine')
 initial={n:sha(ROOT/n) for n in INPUTS}
 results,s=run_python('verify_root15.py',['root15_verified.json']);root=results['root15_verified.json'];steps.append(s)
 results,s=run_python('verify_upper15.py',['upper15_verified.json']);upper=results['upper15_verified.json'];steps.append(s)
 results,s=run_python('anchor_isolation15.py',['anchor_isolation_verified.json']);local=results['anchor_isolation_verified.json'];steps.append(s)
 compiler=shutil.which('g++') or shutil.which('clang++')
 if compiler is None:raise RuntimeError('C++17 compiler and Boost.Multiprecision headers required')
 with (REPORT/'compile.log').open('wb') as f:
  for name in ['enumerate15','audit15']:subprocess.run([compiler,'-O3','-std=c++17',str(ENUM/(name+'.cpp')),'-o',str(RUNTIME/name)],stdout=f,stderr=subprocess.STDOUT,check=True)
 with ThreadPoolExecutor(max_workers=2) as pool:steps.extend(pool.map(family,range(10,15)))
 results,s=run_python('verify_partition15.py',['partition15_verified.json']);part=results['partition15_verified.json'];steps.append(s)
 results,s=run_python('verify_forest15.py',['forest15_verified.json']);forest=results['forest15_verified.json'];steps.append(s)
 results,s=run_python('test_fail_closed15.py',['fail_closed15_verified.json']);negative=results['fail_closed15_verified.json'];steps.append(s)
 rh=initial['core/root15.json'];assert root['root_sha256']==upper['root_sha256']==local['root_sha256']==forest['root_sha256']==part['root_sha256']==rh
 assert forest['local_certificate_sha256']==initial['core/anchor_isolation_certificate.json']
 assert part['residual_sha256']==initial['enumeration/metric15_residuals.json']
 assert root['dimension']==148 and upper['vertices']==26 and upper['edges']==64 and upper['faces']==39 and upper['boundary_caps']==11
 assert upper['exact_equality_center_faces']==13 and len(upper['strict_center_faces'])==4
 assert F(local['anchor_Frobenius_radius'])==F(1,9) and F(local['isolation_margin'])==F(77,540000)>0
 assert part['total']==11950884 and part['residual']==43014 and forest['noncandidate']['cases']==43013
 direct=sum(v['direct'] for v in part['families'].values());farkas=sum(v['farkas'] for v in part['families'].values());excluded=direct+farkas
 assert direct==11729714 and farkas==178156 and excluded==11907870 and excluded+43013+1==part['total']
 assert forest['candidate']['counts']=={'SPLIT':1014,'DUAL':939,'LOCAL':76}
 assert forest['noncandidate']['counts']=={'DUAL':50959,'SPLIT':7946}
 assert forest['candidate']['nodes']==2029 and forest['noncandidate']['nodes']==58905
 assert F(forest['candidate']['minimum_force_margin'])>0 and F(forest['candidate']['minimum_local_margin'])>0 and F(forest['noncandidate']['minimum_force_margin'])>0
 final={n:sha(ROOT/n) for n in INPUTS};assert initial==final,'proof inputs changed during replay'
 out={'verified':True,'theorem':'r_15 = sqrt(t_star), R_15 = 1/sqrt(t_star)','global_optimality_certificate_chain_closed':True,
 'decimal_intervals':root['decimal_intervals'],'kkt_dimension':148,'root_sha256':rh,'covering_upper_verified':True,'enumeration_total':part['total'],
 'metric_direct_excluded':direct,'metric_Farkas_excluded':farkas,'metric_excluded':excluded,'noncandidate_cases':43013,'candidate_cases':1,
 'candidate_tree':forest['candidate'],'noncandidate_forest':forest['noncandidate'],'unresolved_leaves':0,'anchor_isolation_radius':'1/9',
 'labelled_counts_from_independent_recurrence':part['labelled_totals_from_recurrence'],'negative_tests':negative,'input_sha256':initial,'steps':steps,'seconds':time.time()-st,
 'trust_boundary':'Written geometric/topological reductions in PROOF_zh.md; Python integer/Fraction arithmetic; C++17 compiler and Boost integers. All triangulations regenerated, independently BFS-normalized, and checked against rooted counting. No numerical optimizer trusted. Not a Lean/Coq formalization or independent peer review.'}
 temp=MASTER.with_suffix('.tmp');temp.write_text(json.dumps(out,indent=2));os.replace(temp,MASTER)
 print('R15 COMPLETE CERTIFICATE CHAIN VERIFIED',flush=True)
 print(json.dumps({k:out[k] for k in ['verified','theorem','decimal_intervals','enumeration_total','metric_excluded','noncandidate_cases','unresolved_leaves','seconds']},indent=2),flush=True)
 return out
if __name__=='__main__':
 try:main()
 except BaseException:
  MASTER.unlink(missing_ok=True)
  raise
