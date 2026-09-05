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
import hashlib,json,os,shutil,subprocess,sys,time,argparse
OPTIONS=None;CERT=None
INPUTS=['PROOF_zh.md','README.md','verify_all.py',
 *['core/'+n for n in ['system16.py','root16.json','verify_root16.py','interval16.py','exact_arcs.py','verify_upper16.py','anchor_angles16.json','anchor_isolation16.py','anchor_isolation_certificate.json','exact_force16.py','candidate_global16.py','candidate16_tree.json.gz','verify_forest16.py','verify_partition16.py','test_fail_closed16.py']],
 *['enumeration/'+n for n in ['enumerate16.cpp','state_store.hpp','audit16.cpp','metric16_certificate.json','metric16_residuals.json','candidate16_map.json','noncandidate16_trees.jsonl.gz']]]

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
 st=time.time();fam=f'B{B}_I{16-B}';prefix=RUNTIME/fam;out=REPORT/f'audit16_B{B}';log=REPORT/(fam+'.log')
 for p in [out.with_suffix('.json'),out.with_suffix('.residuals')]:p.unlink(missing_ok=True)
 if not OPTIONS.reuse_enumeration:
  prefix.with_suffix('.masks').unlink(missing_ok=True);prefix.with_suffix('.json').unlink(missing_ok=True)
 c=CERT[fam];cs=c['certified_signatures'];lines=[f'9974 21143 100000 1000000000000 {len(cs)}']
 for key,nums in cs.items():
  masks=list(map(int,key.split(','))) if key else [];lines.extend([str(len(masks)),' '.join(map(str,masks)),' '.join(map(str,nums))])
 cuts=prefix.with_suffix('.farkas');cuts.write_text('\n'.join(lines)+'\n')
 with log.open('wb') as f:
  if not OPTIONS.reuse_enumeration:
   subprocess.run([str(RUNTIME/'enumerate16'),str(B),str(16-B),str(prefix)],stdout=f,stderr=subprocess.STDOUT,check=True)
  else:
   f.write(b'REUSE UNTRUSTED ENUMERATION TABLE; ALL ORBITS WILL BE INDEPENDENTLY AUDITED.\n');f.flush()
  assert sha(prefix.with_suffix('.masks'))==c['source_masks_sha256'],'regenerated enumeration differs'
  subprocess.run([str(RUNTIME/'audit16'),str(prefix.with_suffix('.masks')),str(cuts),str(labelled(B,16-B)),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True)
 d=json.load(open(out.with_suffix('.json')));assert d['verified'] is True and int(d['orbit_mass'])==labelled(B,16-B)
 print('PASSED independent BFS census and metric audit',B,16-B,'orbits',d['representatives'],flush=True)
 return {'step':fam,'seconds':time.time()-st,'log_sha256':sha(log),'labelled_count':labelled(B,16-B)}

def main():
 global OPTIONS,CERT
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--reuse-enumeration',action='store_true',help='Reuse untrusted _runtime/*.masks after checking their hashes; still audit every representative independently.')
 parser.add_argument('--jobs',type=int,default=2,help='Concurrent enumeration/audit families (default 2).')
 OPTIONS=parser.parse_args()
 if not 1<=OPTIONS.jobs<=5:raise ValueError('jobs must be between 1 and 5')
 REPORT.mkdir(exist_ok=True);RUNTIME.mkdir(exist_ok=True);MASTER.unlink(missing_ok=True);st=time.time();steps=[]
 if sys.byteorder!='little':raise RuntimeError('The supplied binary interchange format requires a little-endian machine')
 initial={n:sha(ROOT/n) for n in INPUTS}
 CERT=json.load(open(ENUM/'metric16_certificate.json'))
 results,s=run_python('verify_root16.py',['root16_verified.json']);root=results['root16_verified.json'];steps.append(s)
 results,s=run_python('verify_upper16.py',['upper16_verified.json']);upper=results['upper16_verified.json'];steps.append(s)
 results,s=run_python('anchor_isolation16.py',['anchor_isolation_verified.json']);local=results['anchor_isolation_verified.json'];steps.append(s)
 compiler=shutil.which('g++') or shutil.which('clang++')
 if compiler is None:raise RuntimeError('C++17 compiler and Boost.Multiprecision headers required')
 with (REPORT/'compile.log').open('wb') as f:
  for name in ['enumerate16','audit16']:subprocess.run([compiler,'-O3','-std=c++17',str(ENUM/(name+'.cpp')),'-o',str(RUNTIME/name)],stdout=f,stderr=subprocess.STDOUT,check=True)
 with ThreadPoolExecutor(max_workers=OPTIONS.jobs) as pool:steps.extend(pool.map(family,range(11,16)))
 results,s=run_python('verify_partition16.py',['partition16_verified.json']);part=results['partition16_verified.json'];steps.append(s)
 results,s=run_python('verify_forest16.py',['forest16_verified.json']);forest=results['forest16_verified.json'];steps.append(s)
 results,s=run_python('test_fail_closed16.py',['fail_closed16_verified.json']);negative=results['fail_closed16_verified.json'];steps.append(s)
 rh=initial['core/root16.json'];assert root['root_sha256']==upper['root_sha256']==local['root_sha256']==forest['root_sha256']==part['root_sha256']==rh
 assert forest['local_certificate_sha256']==initial['core/anchor_isolation_certificate.json']
 assert part['residual_sha256']==initial['enumeration/metric16_residuals.json']
 assert root['dimension']==140 and upper['vertices']==28 and upper['edges']==69 and upper['faces']==42 and upper['boundary_caps']==12
 assert upper['exact_equality_center_faces']==10 and len(upper['strict_center_faces'])==8
 assert F(local['anchor_Frobenius_radius'])==F(1,10) and F(local['isolation_margin'])==F(43,500000)>0
 assert part['total']==53059205 and part['residual']==176 and forest['noncandidate']['cases']==175
 direct=sum(v['direct'] for v in part['families'].values());farkas=sum(v['farkas'] for v in part['families'].values());excluded=direct+farkas
 assert direct==51464643 and farkas==1594386 and excluded==53059029 and excluded+175+1==part['total']
 assert forest['candidate']['counts']=={'SPLIT':5173,'DUAL':5052,'LOCAL':122}
 assert forest['noncandidate']['counts']=={'DUAL':444,'SPLIT':269}
 assert forest['candidate']['nodes']==10347 and forest['noncandidate']['nodes']==713
 assert F(forest['candidate']['minimum_force_margin'])>0 and F(forest['candidate']['minimum_local_margin'])>0 and F(forest['noncandidate']['minimum_force_margin'])>0
 final={n:sha(ROOT/n) for n in INPUTS};assert initial==final,'proof inputs changed during replay'
 out={'verified':True,'theorem':'r_16 = sqrt(t_star), R_16 = 1/sqrt(t_star)','global_optimality_certificate_chain_closed':True,
 'decimal_intervals':root['decimal_intervals'],'kkt_dimension':140,'root_sha256':rh,'covering_upper_verified':True,'enumeration_total':part['total'],
 'metric_direct_excluded':direct,'metric_Farkas_excluded':farkas,'metric_excluded':excluded,'noncandidate_cases':175,'candidate_cases':1,
 'candidate_tree':forest['candidate'],'noncandidate_forest':forest['noncandidate'],'unresolved_leaves':0,'anchor_isolation_radius':'1/10',
 'labelled_counts_from_independent_recurrence':part['labelled_totals_from_recurrence'],'negative_tests':negative,'input_sha256':initial,'steps':steps,'seconds':time.time()-st,
 'enumeration_mode':('reused_untrusted_tables_independently_audited' if OPTIONS.reuse_enumeration else 'regenerated_and_independently_audited'),
 'external_peer_review_completed':False,'proof_assistant_formalization':False,
 'trust_boundary':'Written geometric/topological reductions in PROOF_zh.md; Python integer/Fraction arithmetic; C++17 compiler and Boost integers. All triangulation tables treated as untrusted, independently BFS-normalized, and checked against rooted counting; default mode regenerates them first. No numerical optimizer trusted. Not a Lean/Coq formalization or independent peer review.'}
 temp=MASTER.with_suffix('.tmp');temp.write_text(json.dumps(out,indent=2));os.replace(temp,MASTER)
 print('R16 COMPLETE CERTIFICATE CHAIN VERIFIED',flush=True)
 print(json.dumps({k:out[k] for k in ['verified','theorem','decimal_intervals','enumeration_total','metric_excluded','noncandidate_cases','unresolved_leaves','seconds']},indent=2),flush=True)
 return out
if __name__=='__main__':
 try:main()
 except BaseException:
  MASTER.unlink(missing_ok=True)
  raise
