#!/usr/bin/env python3
"""Fail-closed replay of the finite certificates in the r_13 proof.
Python standard library + a C++17 compiler with Boost.Multiprecision headers.
Every checker is run anew. Existing result JSON files are never accepted alone.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,sys,json,time,subprocess,shutil,gzip,hashlib
if not __debug__ or sys.flags.optimize:
 raise RuntimeError('Run with ordinary Python, never -O/-OO/PYTHONOPTIMIZE.')
ROOT=Path(__file__).resolve().parent;CORE=ROOT/'core';ENUM=ROOT/'enumeration';WORK=ROOT/'work';REPORT=ROOT/'reports'
WORK.mkdir(exist_ok=True);REPORT.mkdir(exist_ok=True)
ENV=dict(os.environ,PYTHONOPTIMIZE='0',PYTHONDONTWRITEBYTECODE='1')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name,args):
 st=time.time();p=REPORT/(name+'.log')
 with p.open('w',encoding='utf8') as log:
  r=subprocess.run(args,cwd=ROOT,env=ENV,stdout=log,stderr=subprocess.STDOUT)
 if r.returncode:
  print(p.read_text()[-8000:]);raise RuntimeError(f'{name}: return code {r.returncode}')
 print(name,'PASS',flush=True);return {'check':name,'returncode':0,'seconds':time.time()-st,'log_sha256':sha(p)}
def read_verified(p):
 x=json.loads(p.read_text());assert x.get('verified') is True;return x

def main():
 st=time.time();out=REPORT/'MASTER_VERIFIED.json'
 # A failed replay must not leave an old successful master result at this location.
 if out.exists():out.unlink()
 records=[]
 manifest=ROOT/'MANIFEST_SHA256.json'
 if manifest.exists():
  mm=json.loads(manifest.read_text())
  for rel,h in mm.items():
   p=ROOT/rel;assert p.is_file() and sha(p)==h,('manifest mismatch',rel)
  print('input_manifest PASS',flush=True)
 for name in ['verify_interfaces','verify_root_repaired','verify_upper','anchor_isolation','candidate_global','verify_partition','noncandidate_global','test_fail_closed']:
  records.append(run(name,[sys.executable,'-B',str(CORE/(name+'.py'))]))
 it=read_verified(CORE/'interfaces_verified.json');rt=read_verified(CORE/'cover13_krawczyk_core_verified.json')
 up=read_verified(CORE/'cover13_upper_verified.json');loc=read_verified(CORE/'anchor_isolation_verified.json')
 cg=read_verified(CORE/'candidate_global_verified.json');pt=read_verified(CORE/'partition_verified.json');ng=read_verified(ENUM/'noncandidate_verified.json')
 assert rt['dimension']==it['dimension']==119
 assert rt['t_interval']==up['t_interval']==it['t_root_interval']
 assert up['krawczyk_cert_sha256']==loc['root_sha256']==sha(CORE/'cover13_krawczyk_cert.json')
 assert pt['residual']==ng['cases']+1==1761
 assert cg['leaves']['DUAL']+cg['leaves']['LOCAL']==916 and cg['nodes']==1831
 assert ng['leaf_counts']=={'DUAL':2334} and ng['total_nodes']==2908
 compiler=shutil.which(os.environ.get('CXX','g++')) or shutil.which('clang++')
 if not compiler:raise RuntimeError('A C++17 compiler (g++ or clang++) is required.')
 exe=WORK/('orbit_audit.exe' if os.name=='nt' else 'orbit_audit')
 records.append(run('compile_orbit_audit',[compiler,'-O3','-std=c++17',str(ENUM/'audit.cpp'),'-o',str(exe)]))
 cert=json.loads((ENUM/'metric_partition_certificate.json').read_text())
 def family(B):
  I=13-B;fam=f'B{B}_I{I}';path=WORK/(fam+'.txt');path.write_bytes(gzip.decompress((ENUM/(fam+'.txt.gz')).read_bytes()))
  assert sha(path)==cert[fam]['source_txt_sha256']
  rec=run('orbit_'+fam,[str(exe),str(path),str(REPORT/(fam+'_orbit_verified.json'))])
  v=read_verified(REPORT/(fam+'_orbit_verified.json'));assert (v['B'],v['I'])==(B,I)
  assert v['orbit_mass']==it['labeled_counts_by_recurrence'][fam]
  assert v['representatives']==sum(pt['families'][fam].values())
  return rec,v
 # The audits are independent. Running them in parallel affects no arithmetic.
 with ThreadPoolExecutor(max_workers=4) as pool:
  results=list(pool.map(family,range(9,13)))
 orbits={}
 for rec,v in results:records.append(rec);orbits[f"B{v['B']}_I{v['I']}"]=v
 total=sum(v['representatives'] for v in orbits.values());assert total==pt['total']==308198
 excluded=sum(v.get('direct',0)+v.get('farkas',0) for v in pt['families'].values());assert excluded+ng['cases']+1==total
 report={'verified':True,'finite_certificates_verified':True,'enumerated_cases_complete':True,
 'theorem':'r_13 = sqrt(t_star), with t_star defined by the enclosed 119-variable polynomial root certificate.',
 'mathematical_reduction':'See PROOF_zh.md: irreducible-cover no-return lemma, dual triangulation, stellar padding, root-edge counting, and exact dual/local lemmas.',
 'r_strict_interval':it['r_strict_interval'],'R_strict_interval':it['R_strict_interval'],
 'total_topology_orbits':total,'metric_excluded':excluded,'noncandidate_residuals_excluded':ng['cases'],'candidate_orbits':1,
 'candidate_tree_nodes':cg['nodes'],'candidate_force_leaves':cg['leaves']['DUAL'],'candidate_local_leaves':cg['leaves']['LOCAL'],
 'noncandidate_tree_nodes':ng['total_nodes'],'noncandidate_force_leaves':ng['leaf_counts']['DUAL'],
 'candidate_minimum_radius_margin':cg['minimum_dual_radius_margin'],'noncandidate_minimum_radius_margin':ng['minimum_radius_margin'],
 'orbit_audits':orbits,'checks':records,'seconds':time.time()-st,
 'status_note':'This is an exact-arithmetic computer-assisted proof with written geometric reductions, not a Lean/Coq formalization or a claim of peer review.'}
 out.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
