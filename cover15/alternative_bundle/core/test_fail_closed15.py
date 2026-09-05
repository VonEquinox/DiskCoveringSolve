"""Adversarial regression tests supplementing full exact replay."""
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import copy,gzip,json,shutil,struct,subprocess,sys,tempfile,time
import exact_force15 as G
from verify_forest15 import traverse
R=Path(__file__).resolve().parent;P=R.parent;E=P/'enumeration'

def rejected(fn):
 try:fn()
 except (AssertionError,ValueError,KeyError,IndexError,RuntimeError,StopIteration):return True
 return False

def main():
 st=time.time();tests={}
 with gzip.open(E/'noncandidate15_trees.jsonl.gz','rt',encoding='utf8') as f:header=json.loads(f.readline());case=json.loads(f.readline())
 B=case['B'];faces=case['faces'];assert len(case['nodes'])==1 and case['nodes'][0]['kind']=='DUAL';traverse(B,faces,case['nodes'],False)
 bad=copy.deepcopy(case['nodes']);bad[0]['force_num'][0][0]+=1
 assert rejected(lambda:traverse(B,faces,bad,False));tests['one_integer_force_unit_corrupted']='rejected'
 assert rejected(lambda:traverse(B,faces,[{'kind':'UNRESOLVED'}],False));tests['unresolved_leaf']='rejected'
 bad=copy.deepcopy(case['nodes'])+[{'kind':'EMPTY'}]
 assert rejected(lambda:traverse(B,faces,bad,False));tests['unreachable_extra_node']='rejected'
 M=G.bounds(B,faces);lo,hi=G.rootbox(M);lo,hi=G.tighten(M,lo,hi);mid=(lo[0]+hi[0])/2
 bad=[{'kind':'SPLIT','axis':0,'mid':str(mid),'children':[1,1]},{'kind':'EMPTY'}]
 assert rejected(lambda:traverse(B,faces,bad,False));tests['duplicate_split_child']='rejected'
 assert rejected(lambda:traverse(B,faces,[{'kind':'LOCAL'}],False));tests['local_leaf_in_noncandidate_case']='rejected'
 with tempfile.TemporaryDirectory(prefix='r15-reject-') as td:
  T=Path(td)
  for n in ['system15.py','verify_root15.py','root15.json']:shutil.copyfile(R/n,T/n)
  c=json.load(open(T/'root15.json'));c['xnum'][0]=str(int(c['xnum'][0])+int(c['Qx'])//1000);(T/'root15.json').write_text(json.dumps(c))
  rr=subprocess.run([sys.executable,'-S','-B',str(T/'verify_root15.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0;tests['corrupted_algebraic_root']='rejected'
  with (P/'_runtime/B10_I5.masks').open('rb') as f:head=f.read(20);hh=list(struct.unpack('<5I',head));rec=f.read(hh[3]*2)
  hh[4]=2;(T/'duplicate.masks').write_bytes(struct.pack('<5I',*hh)+rec+rec)
  (T/'empty.farkas').write_text('10306 21954 100000 1000000000000 0\n')
  rr=subprocess.run([str(P/'_runtime/audit15'),str(T/'duplicate.masks'),str(T/'empty.farkas'),'15274459200',str(T/'duplicate_audit')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'duplicate orbit' in rr.stdout;tests['duplicate_orbit_representative']='rejected'
  hh[4]=1;(T/'invalid.masks').write_bytes(struct.pack('<5I',*hh)+b'\x00\x00'+rec[2:])
  rr=subprocess.run([str(P/'_runtime/audit15'),str(T/'invalid.masks'),str(T/'empty.farkas'),'15274459200',str(T/'invalid_audit')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'bad triangle' in rr.stdout;tests['invalid_triangle_mask']='rejected'
  # Keep a correct global header but remove all but one noncandidate tree.
  C=T/'missing_case'/'core';EN=T/'missing_case'/'enumeration';C.mkdir(parents=True);EN.mkdir()
  for n in ['system15.py','root15.json','interval15.py','exact_arcs.py','exact_force15.py','candidate_global15.py','verify_forest15.py','candidate15_tree.json.gz','anchor_isolation_certificate.json']:shutil.copyfile(R/n,C/n)
  for n in ['metric15_residuals.json','candidate15_map.json']:shutil.copyfile(E/n,EN/n)
  with gzip.open(EN/'noncandidate15_trees.jsonl.gz','wt',encoding='utf8') as f:f.write(json.dumps(header)+'\n'+json.dumps(case)+'\n')
  rr=subprocess.run([sys.executable,'-S','-B',str(C/'verify_forest15.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'AssertionError' in rr.stdout;tests['missing_noncandidate_cases']='rejected'
 rr=subprocess.run([sys.executable,'-S','-B','-O',str(P/'verify_all.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 assert rr.returncode!=0 and b'Run without' in rr.stdout;tests['assertions_disabled']='rejected'
 out={'verified':True,'tests':tests,'seconds':time.time()-st};(R/'fail_closed15_verified.json').write_text(json.dumps(out,indent=2));print('NEGATIVE TESTS PASSED',tests,flush=True);return out
if __name__=='__main__':main()
