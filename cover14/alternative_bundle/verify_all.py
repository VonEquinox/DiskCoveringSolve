#!/usr/bin/env python3
"""Full, fail-closed n=14 proof replay. Requires only Python stdlib, C++17 and Boost headers.
No previous VERIFIED file is an input. Every successful run rechecks every gate.
"""
if not __debug__:
 raise RuntimeError('Exact verification requires Python assertions; do not use -O/-OO/PYTHONOPTIMIZE')
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import sys,json,subprocess,time,gzip,tempfile,hashlib,shutil,os
BASE=Path(__file__).resolve().parent
LOG=BASE/'replay_logs';LOG.mkdir(exist_ok=True)
REPORT=BASE/'r14_MASTER_VERIFIED.json'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def decimal(s):
 v=F(s)
 with localcontext() as ctx:
  ctx.prec=110;return str(Decimal(v.numerator)/Decimal(v.denominator))

def main():
 start=time.time();stages={};audits={};out={'verified':False,'stages':stages}
 REPORT.write_text(json.dumps(out,indent=2))
 def stage(name,args,result):
  t0=time.time();p=BASE/result
  if p.exists():p.unlink()
  with open(LOG/(name+'.log'),'w') as f:r=subprocess.run(args,cwd=BASE,stdout=f,stderr=subprocess.STDOUT)
  assert r.returncode==0,(name,'process failed',r.returncode)
  dat=json.loads(p.read_text());assert dat.get('verified') is True,(name,'result not verified')
  stages[name]={'verified':True,'seconds':time.time()-t0,'report':result};print(name,'VERIFIED',flush=True);return dat
 try:
  root=stage('root',[sys.executable,'-B','verify_root.py'],'root_verified.json')
  sy=stage('symmetry',[sys.executable,'-B','verify_symmetry.py'],'symmetry_verified.json')
  up=stage('upper',[sys.executable,'-B','verify_upper.py'],'upper_verified.json')
  loc=stage('anchor_isolation',[sys.executable,'-B','anchor_isolation.py'],'anchor_isolation_verified.json')
  part=stage('metric_partition',[sys.executable,'-B','verify_partition.py'],'partition_verified.json')
  # A fresh executable, and precisely the decompressed files used by the metric verifier.
  from verify_partition import triangulation_count
  with tempfile.TemporaryDirectory(prefix='r14_exact_audit_') as tmp:
   tmp=Path(tmp);exe=tmp/'audit'
   cmd=['g++','-O3','-std=c++17',str(BASE/'audit.cpp'),'-o',str(exe)]
   with open(LOG/'compile.log','w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
   assert r.returncode==0,('C++17/Boost compiler failed',r.returncode)
   for B in range(10,14):
    I=14-B;name=f'B{B}_I{I}';raw=gzip.decompress((BASE/'enumeration'/(name+'.txt.gz')).read_bytes());src=tmp/(name+'.txt');src.write_bytes(raw)
    assert hashlib.sha256(raw).hexdigest()==part['families'][name]['source_txt_sha256']
    a=stage('orbit_audit_'+name,[str(exe),str(src),str(BASE/'enumeration'/('audit_'+name+'.json'))],'enumeration/audit_'+name+'.json')
    assert a['B']==B and a['I']==I and a['representatives']==part['families'][name]['orbits']
    assert int(a['orbit_mass'])==triangulation_count(B,I)==part['families'][name]['labeled_count_by_root_peeling']
    audits[name]=a
  non=stage('all_noncandidates',[sys.executable,'-B','verify_force.py'],'noncandidate_verified.json')
  cand=stage('candidate_domain',[sys.executable,'-B','verify_force.py','candidate'],'candidate_verified.json')
  neg=stage('negative_tests',[sys.executable,'-B','test_fail_closed.py'],'negative_tests_verified.json')
  assert part['candidate_residual_index']==cand['residual_index']==685
  assert non['cases']+1==part['residuals'] and part['metric_excluded']+non['cases']+1==part['total_orbits']
  assert non['residual_sha256']==part['residual_sha256']
  h=digest(BASE/'root_certificate.json')
  assert root['root_sha256']==sy['root_sha256']==up['root_sha256']==loc['root_sha256']==cand['root_sha256']==h
  assert all(v['verified'] is True for v in stages.values())
  readonly=list(BASE.glob('*.py'))+[BASE/'audit.cpp',BASE/'root_certificate.json',BASE/'anchor_isolation_certificate.json',BASE/'candidate_certificate.json.gz',BASE/'noncandidate_certificate.json.gz']+list((BASE/'enumeration').glob('*.txt.gz'))+[BASE/'enumeration/metric_partition_certificate.json',BASE/'enumeration/metric_residuals.json']
  out={'verified':True,'theorem':'r_14 = sqrt(t_star), defined by the isolated 116-variable KKT system and rational box in root_certificate.json','r_strict_interval':[decimal(v) for v in root['r_strict_interval']],'R_strict_interval':[decimal(v) for v in root['R_strict_interval']],'stages':stages,'enumeration':audits,'partition':{'all_orbits':part['total_orbits'],'metric_excluded':part['metric_excluded'],'noncandidate_cases':non['cases'],'candidate_cases':1,'unresolved_cases':0},'noncandidate_trees':{'nodes':non['nodes'],'counts':non['counts'],'min_radius_margin':non['minimum_radius_margin']},'candidate_tree':{'nodes':cand['nodes'],'counts':cand['counts'],'max_local_anchor_distance_squared':cand['max_local_anchor_distance_squared']},'local_theorem':{'gamma':loc['gamma'],'lambda':loc['lambda'],'nu':loc['nu'],'anchor_radius':loc['anchor_radius'],'isolation_margin':loc['isolation_margin']},'upper_structure':{'center_faces':up['center_faces'],'positive_stress_faces':up['positive_stress_faces'],'extra_tight_faces':up['extra_tight_faces'],'strict_slack_faces':up['strict_slack_faces'],'omitted_degenerate_radial_faces':up['omitted_degenerate_radial_faces']},'negative_tests':neg['tests_rejected'],'input_sha256':{str(p.relative_to(BASE)):digest(p) for p in sorted(readonly)},'seconds':time.time()-start,'scope':'Written geometric/combinatorial reductions plus exact executable finite certificates. Not third-party peer review or a Lean/Coq formalization.'}
  REPORT.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True);return 0
 except Exception as e:
  out={'verified':False,'stages':stages,'error':repr(e),'seconds':time.time()-start};REPORT.write_text(json.dumps(out,indent=2));raise
if __name__=='__main__':raise SystemExit(main())
