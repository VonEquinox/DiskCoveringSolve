#!/usr/bin/env python3
"""Small negative tests: proof inputs must not silently pass when incomplete/corrupt."""
from pathlib import Path
from copy import deepcopy
from tempfile import TemporaryDirectory
import gzip,json,shutil,subprocess,sys
if not __debug__:raise RuntimeError('Do not disable assertions')
import noncandidate_global as N
import candidate_global as C
ROOT=Path(__file__).resolve().parent;ENUM=ROOT.parent/'enumeration'

def must_reject(call,label):
 try:call()
 except (AssertionError,ValueError,KeyError,FileNotFoundError):return label
 raise RuntimeError(label+' was unexpectedly accepted')

def main():
 tests=[]
 with gzip.open(ENUM/'noncandidate_certificate.json.gz','rt') as f:cert=json.load(f)
 first=cert['cases'][0];assert len(first['nodes'])==1 and first['nodes'][0]['kind']=='DUAL'
 B=first['B'];ed=N.edges(B,first['faces']);M=N.bounds(B,first['faces']);lo,hi=N.tighten(M,*N.rootbox(M))
 assert N.check_dual(B,ed,lo,hi,first['nodes'][0])>0
 bad=deepcopy(first['nodes'][0]);bad['force_num'][0][0]+=1
 tests.append(must_reject(lambda:N.check_dual(B,ed,lo,hi,bad),'one-unit force corruption rejected'))
 with TemporaryDirectory() as td:
  p=Path(td)
  for name in ['metric_residuals_new.json','candidate_orbit_map.json']:shutil.copy2(ENUM/name,p/name)
  empty=dict(cert,cases=[])
  with gzip.open(p/'noncandidate_certificate.json.gz','wt') as f:json.dump(empty,f)
  old=N.OUT;N.OUT=p
  try:tests.append(must_reject(N.verify,'empty noncandidate certificate set rejected'))
  finally:N.OUT=old
 with TemporaryDirectory() as td:
  p=Path(td)
  for name in ['cover13_krawczyk_cert.json','cover13_kkt_root_119d.json','anchor_isolation_certificate.json']:shutil.copy2(ROOT/name,p/name)
  with gzip.open(ROOT/'candidate_global_certificate.json.gz','rt') as f:cc=json.load(f)
  cc['nodes'][0]={'kind':'UNRESOLVED'}
  with gzip.open(p/'candidate_global_certificate.json.gz','wt') as f:json.dump(cc,f)
  old=C.ROOT;C.ROOT=p
  try:tests.append(must_reject(C.verify,'unresolved candidate leaf rejected'))
  finally:C.ROOT=old
 z=subprocess.run([sys.executable,'-O',str(ROOT.parent/'verify_all.py')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 assert z.returncode!=0 and b'Run with ordinary Python' in z.stderr
 tests.append('optimized-Python bypass rejected')
 out={'verified':True,'tests':tests};(ROOT/'negative_tests_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
