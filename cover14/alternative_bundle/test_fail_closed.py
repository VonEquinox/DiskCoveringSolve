#!/usr/bin/env python3
"""Adversarial rejection tests of the accepting checkers, not checks of cached success flags."""
from pathlib import Path
import json,gzip,copy,subprocess,sys,time
from verify_force import *
if not __debug__:raise RuntimeError('Assertions must be enabled')
BASE=Path(__file__).resolve().parent

def main():
 st=time.time();passed=[]
 def reject(name,fn):
  try:fn()
  except (AssertionError,RuntimeError,KeyError,ValueError,TypeError,IndexError):passed.append(name);return
  raise AssertionError(('bad input was accepted',name))
 with gzip.open(BASE/'candidate_certificate.json.gz','rt') as f:ca=json.load(f)
 with gzip.open(BASE/'noncandidate_certificate.json.gz','rt') as f:nc=json.load(f)
 one=next(c for c in nc['cases'] if c['nodes'][0]['kind']=='DUAL');B=one['B'];ed=edges(B,one['faces']);M=bounds(B,one['faces']);lo,hi=rootbox(M);lo,hi=tighten(M,lo,hi)
 rec=copy.deepcopy(one['nodes'][0]);rec['force_num'][0][0]+=1
 reject('one-unit force corruption',lambda:check_dual(B,ed,lo,hi,rec))
 rec2=copy.deepcopy(one['nodes'][0]);rec2['lambda_num'][0]=-1
 reject('negative halfplane multiplier',lambda:check_dual(B,ed,lo,hi,rec2))
 reject('empty noncandidate partition',lambda:check_residual_partition([],7232))
 keys=[c['residual_index'] for c in nc['cases']];reject('one missing residual case',lambda:check_residual_partition(keys[1:],7232))
 bad=copy.deepcopy(ca);bad['nodes'][0]={'kind':'PENDING'};reject('unresolved candidate root',lambda:replay_tree(bad,True))
 bad2=copy.deepcopy(ca);bad2['nodes'].append({'kind':'EMPTY'});reject('unreachable extra leaf',lambda:replay_tree(bad2,True))
 bad3=copy.deepcopy(ca);bad3['nodes'][0]['children']=[1,1];reject('duplicate split children',lambda:replay_tree(bad3,True))
 r=subprocess.run([sys.executable,'-O',str(BASE/'verify_all.py')],capture_output=True,text=True)
 assert r.returncode!=0 and 'assert' in (r.stdout+r.stderr).lower();passed.append('disabled Python assertions')
 out={'verified':True,'tests_rejected':passed,'number':len(passed),'seconds':time.time()-st};(BASE/'negative_tests_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
