"""Small adversarial regressions; these supplement, not replace, proof replay."""
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import copy,gzip,json,shutil,subprocess,sys,tempfile,time
import exact_force14 as G
from verify_forest14 import traverse
R=Path(__file__).resolve().parent;P=R.parent;E=P/'enumeration'

def rejected(fn):
 try:fn()
 except (AssertionError,ValueError,KeyError,IndexError,RuntimeError):return True
 return False

def main():
 st=time.time();tests={}
 with gzip.open(E/'noncandidate14_trees.json.gz','rt') as f:data=json.load(f)
 case=data['cases'][0];B=case['B'];faces=case['faces'];assert len(case['nodes'])==1 and case['nodes'][0]['kind']=='DUAL'
 traverse(B,faces,case['nodes'],False)
 bad=copy.deepcopy(case['nodes']);bad[0]['force_num'][0][0]+=1
 assert rejected(lambda:traverse(B,faces,bad,False));tests['one_integer_force_unit_corrupted']='rejected'
 assert rejected(lambda:traverse(B,faces,[{'kind':'UNRESOLVED'}],False));tests['unresolved_leaf']='rejected'
 bad=copy.deepcopy(case['nodes'])+[{'kind':'EMPTY'}]
 assert rejected(lambda:traverse(B,faces,bad,False));tests['unreachable_extra_node']='rejected'
 M=G.bounds(B,faces);lo,hi=G.rootbox(M);lo,hi=G.tighten(M,lo,hi);mid=(lo[0]+hi[0])/2
 bad=[{'kind':'SPLIT','axis':0,'mid':str(mid),'children':[1,1]},{'kind':'EMPTY'}]
 assert rejected(lambda:traverse(B,faces,bad,False));tests['duplicate_split_child']='rejected'
 with tempfile.TemporaryDirectory(prefix='r14-reject-') as td:
  T=Path(td)
  for n in ['system14.py','verify_root14.py','root14.json']:shutil.copyfile(R/n,T/n)
  c=json.load(open(T/'root14.json'));c['xnum'][0]=str(int(c['xnum'][0])+int(c['Qx'])//1000);(T/'root14.json').write_text(json.dumps(c))
  rr=subprocess.run([sys.executable,'-S','-B',str(T/'verify_root14.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0;tests['corrupted_algebraic_root']='rejected'
  with gzip.open(E/'B10_I4.txt.gz','rt') as f:head=f.readline();line=f.readline()
  (T/'duplicate.txt').write_text('10 4 16 2\n'+line+line)
  rr=subprocess.run([str(P/'_runtime/audit14'),str(T/'duplicate.txt')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0;tests['duplicate_orbit_representative']='rejected'
 rr=subprocess.run([sys.executable,'-S','-B','-O',str(P/'verify_all.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 assert rr.returncode!=0;tests['assertions_disabled']='rejected'
 out={'verified':True,'tests':tests,'seconds':time.time()-st};(R/'fail_closed14_verified.json').write_text(json.dumps(out,indent=2));print('NEGATIVE TESTS PASSED',tests,flush=True)
 return out
if __name__=='__main__':main()
