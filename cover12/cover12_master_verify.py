#!/usr/bin/env python3
"""Master replay for the computer-assisted proof of the n=12 disk cover.

Every mathematical acceptance decision is delegated to a verifier using exact
integer/Fraction arithmetic, explicit Taylor bounds, interval arithmetic, or
symbolic polynomial identities. Floating-point discovery files are never
accepted without re-rationalization and exact checking.
"""
from __future__ import annotations
from pathlib import Path
import subprocess,sys,time,json,hashlib,os
ROOT=Path(__file__).resolve().parent
PY=sys.executable

MODULES=[
 ('symmetric_kkt',[PY,str(ROOT/'cover12_symmetric_kkt_exact.py'),'verify',str(ROOT/'cover12_symmetric_kkt_cert.json')],ROOT/'cover12_symmetric_kkt_verified.json'),
 ('kkt_structure',[PY,str(ROOT/'cover12_kkt_structure_audit.py')],ROOT/'cover12_kkt_structure_audit_verified.json'),
 ('thresholds',[PY,str(ROOT/'cover12_threshold_verify.py')],ROOT/'cover12_threshold_verified.json'),
 ('upper_cover',[PY,str(ROOT/'cover12_upper_exact.py')],ROOT/'cover12_upper_exact_verified.json'),
 ('metric_farkas',[PY,str(ROOT/'metric_exact.py'),'verify',str(ROOT/'metric_exact_cert_safe.json')],ROOT/'metric_exact_safe_verified.json'),
 ('topology_partition',[PY,str(ROOT/'cover12_topology_partition_verify.py')],ROOT/'cover12_topology_partition_verified.json'),
 ('residual_energies',[PY,str(ROOT/'residual_energy_exact_safe.py'),'verify',str(ROOT/'residual_energy_cert_safe.json')],ROOT/'residual_energy_safe_verified.json'),
 ('candidate_convexity',[PY,str(ROOT/'cover12_candidate_convex_exact.py')],ROOT/'cover12_candidate_convex_verified.json'),
]

def sha256(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 t0=time.time();runs=[];logdir=ROOT/'master_logs';logdir.mkdir(exist_ok=True)
 env=dict(os.environ);env['PYTHONHASHSEED']='0'
 for name,cmd,report in MODULES:
  st=time.time();log=logdir/f'{name}.log'
  with log.open('w') as f:
   p=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,env=env)
  if p.returncode:
   print(log.read_text()[-8000:]);raise SystemExit(f'{name} failed with code {p.returncode}')
  d=json.load(report.open());assert d.get('verified') is True,(name,d)
  runs.append({'name':name,'seconds':time.time()-st,'report':report.name,'log':str(log.relative_to(ROOT))})
  print('VERIFIED',name,'seconds',runs[-1]['seconds'],flush=True)
 # Cross-module invariants, intentionally repeated at the master level.
 k=json.load((ROOT/'cover12_symmetric_kkt_verified.json').open())
 th=json.load((ROOT/'cover12_threshold_verified.json').open())
 up=json.load((ROOT/'cover12_upper_exact_verified.json').open())
 met=json.load((ROOT/'metric_exact_safe_verified.json').open())
 part=json.load((ROOT/'cover12_topology_partition_verified.json').open())
 res=json.load((ROOT/'residual_energy_safe_verified.json').open())
 cand=json.load((ROOT/'cover12_candidate_convex_verified.json').open())
 struct=json.load((ROOT/'cover12_kkt_structure_audit_verified.json').open())
 assert k['contraction_inf_norm']<1
 assert met['total_orbits']==38208 and met['total_residual']==132
 assert sum(x['labeled'] for x in met['cases'])==2831335
 assert part['safe_residuals']==132 and part['noncandidate_exact']==131
 assert res['topologies']==131 and res['partition_verified'] and res['min_margin']>0
 assert cand['candidate_orbit']==21015 and cand['min_ldl_pivot']>0
 assert up['combinatorial_disk_verified'] and up['faces']==31 and up['min_inactive_r2_margin']>0
 assert struct['active_edges']==36 and struct['force_components']==36 and struct['anchor_tangential_components']==9
 assert th['T_minus_t_upper_float']>0 and th['sin2_delta_half_minus_T_float']>0
 # Decisive input hashes. Reports/logs are omitted because they are regenerated.
 inputs=[
  'cover12_symmetric_kkt_exact.py','cover12_symmetric_kkt_cert.json',
  'cover12_kkt_structure_audit.py','cover12_threshold_verify.py','cover12_upper_exact.py','cover12_topology_partition_verify.py',
  'metric_exact.py','metric_exact_cert_safe.json','residual_energy_exact_safe.py',
  'residual_energy_cert_safe.json','cover12_candidate_convex_exact.py','fixed_trig.py','exact_core.py',
  'enumerate_orbits.py',
 ]+[f'triangulations_B{B}_I{I}_orbits.pkl' for B,I in [(9,1),(9,2),(9,3),(10,1),(10,2),(11,1)]]
 hashes={p:sha256(ROOT/p) for p in inputs}
 out={
  'verified':True,
  'theorem':'r_12 equals the r-coordinate of the unique D3 candidate root isolated by cover12_symmetric_kkt_cert.json',
  'r_interval':k['r_interval'],'t_interval':k['t_interval'],
  'decimal_r_center':k['center_decimal']['r'],
  'r_decimal_high_precision':'0.361102963744508644113087702017065180848530565591618515449660087986314559381989395...',
  'T':th['T'],'T_minus_t_upper':th['T_minus_t_upper'],
  'orbits':met['total_orbits'],'labeled_triangulations':sum(x['labeled'] for x in met['cases']),
  'farkas_excluded':met['total_orbits']-met['total_residual'],'residual_topologies':met['total_residual'],
  'noncandidate_exact':res['topologies'],'candidate_orbit':[9,3,21015],
  'minimum_noncandidate_margin':res['min_margin'],'candidate_min_ldl_pivot':cand['min_ldl_pivot'],
  'upper_inactive_margin':up['min_inactive_r2_margin'],
  'modules':runs,'input_sha256':hashes,'seconds':time.time()-t0,
 }
 (ROOT/'cover12_master_verified.json').write_text(json.dumps(out,indent=2))
 print(json.dumps({k:v for k,v in out.items() if k not in ('input_sha256','r_interval','t_interval','T_minus_t_upper')},indent=2))
if __name__=='__main__':main()
