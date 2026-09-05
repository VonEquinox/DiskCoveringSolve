"""Deliberate corruptions must be rejected by actual accepting functions."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import contextlib,copy,gzip,io,json,shutil,subprocess,sys,tempfile,time
import exact_force15 as G
import verify_forest15 as V
import anchor_isolation15 as A
import system15 as S
R=Path(__file__).resolve().parent;ROOT=R.parent;ENUM=ROOT/'enumeration';WORK=ROOT/'replay_workspace'

def main():
 st=time.time();tests=[]
 def reject(name,fn):
  try:fn()
  except (AssertionError,ValueError,RuntimeError,KeyError,IndexError,TypeError):tests.append(name);return
  raise AssertionError('Corruption was not rejected: '+name)
 def nonzero(name,cmd,cwd=None,expected=None):
  p=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  assert p.returncode!=0,(name,p.stdout)
  if expected:assert expected in p.stdout,(name,p.stdout)
  tests.append(name)
 with gzip.open(R/'candidate15_tree.json.gz','rt') as f:t=json.load(f)
 B=t['B'];faces=t['faces'];nodes=t['nodes']
 reject('unknown unresolved leaf',lambda:V.traverse(B,faces,[{'kind':'UNRESOLVED'}],True))
 reject('whole domain falsely marked local',lambda:V.traverse(B,faces,[{'kind':'LOCAL'}],True))
 bad=copy.deepcopy(nodes);assert bad[0]['kind']=='SPLIT';bad[0]['children'][1]=bad[0]['children'][0]
 reject('duplicate branch child',lambda:V.traverse(B,faces,bad,True))
 bad=copy.deepcopy(nodes);bad[0]['children'][0]=0
 reject('cyclic branch edge',lambda:V.traverse(B,faces,bad,True))
 bad=copy.deepcopy(nodes);bad.append({'kind':'EMPTY'})
 reject('unreachable appended node',lambda:V.traverse(B,faces,bad,True))
 # Find a genuine accepted force leaf, retaining its exact propagated box.
 M=G.bounds(B,faces);lo,hi=G.rootbox(M);stack=[(0,lo,hi)]
 while stack:
  k,l,h=stack.pop();l,h=G.tighten(M,l,h);n=nodes[k]
  if n['kind']=='DUAL':break
  if n['kind']=='SPLIT':
   j=n['axis'];m=F(n['mid']);lh=list(h);lh[j]=m;rl=list(l);rl[j]=m
   stack.extend([(n['children'][1],rl,list(h)),(n['children'][0],list(l),lh)])
 else:raise AssertionError('No force leaf')
 ed=G.edges(B,faces);assert G.check(B,ed,l,h,n)>0
 bad=copy.deepcopy(n);bad['force_num'][0][0]+=1
 reject('one integer force unit changed',lambda:G.check(B,ed,l,h,bad))
 bad=copy.deepcopy(n);bad['lambda_num'][0]=-1
 reject('negative halfplane multiplier',lambda:G.check(B,ed,l,h,bad))
 bad=copy.deepcopy(n);bad['force_num']=bad['force_num'][:-1]
 reject('missing rod force',lambda:G.check(B,ed,l,h,bad))
 rows=json.load(open(ENUM/'metric15_residuals.json'));cm=json.load(open(ENUM/'candidate15_map.json'));cand=V.check_mapping(rows)
 allidx=set(range(len(rows)))-{cand};missing=allidx-{min(allidx)}
 reject('missing residual topology',lambda:V.check_complete_indices(missing,rows,cand))
 reject('candidate duplicated in residual forest',lambda:V.check_complete_indices(allidx|{cand},rows,cand))
 with tempfile.TemporaryDirectory(prefix='r15_reject_') as td:
  q=Path(td)
  wrong=copy.deepcopy(cm);wrong['map'][0],wrong['map'][1]=wrong['map'][1],wrong['map'][0]
  (q/'candidate15_map.json').write_text(json.dumps(wrong))
  old=V.ENUM
  try:
   V.ENUM=q;reject('incorrect candidate graph permutation',lambda:V.check_mapping(rows))
  finally:V.ENUM=old
  shutil.copy2(R/'root15.json',q/'root15.json')
  shutil.copy2(R/'system15.py',q/'system15.py')
  mat=json.load(open(R/'anchor_isolation_certificate.json'));mat['Rnum'][0][0]=0
  (q/'anchor_isolation_certificate.json').write_text(json.dumps(mat))
  old=A.ROOT
  try:
   A.ROOT=q
   with contextlib.redirect_stdout(io.StringIO()):reject('singular alleged PD congruence',A.main)
  finally:A.ROOT=old
  shutil.copy2(R/'system15.py',q/'system15.py');shutil.copy2(R/'verify_root15.py',q/'verify_root15.py')
  root=json.load(open(q/'root15.json'));root['xnum'][0]=str(int(root['xnum'][0])+int(root['Qx'])//10)
  (q/'root15.json').write_text(json.dumps(root))
  nonzero('algebraic root midpoint displaced',[sys.executable,'-S','-B',q/'verify_root15.py'],q,'AssertionError')
  with open(WORK/'B10_I5.txt') as f:
   header=f.readline();first=f.readline()
  (q/'duplicate.txt').write_text('10 5 18 2\n'+first+first)
  nonzero('duplicate valid enumeration orbit',[WORK/'audit15',q/'duplicate.txt',q/'bad_audit.json'],expected='duplicate orbit')
  (q/'one.txt').write_text('10 5 18 1\n'+first)
  (q/'negative.classes').write_text('10 1\n1 0 -1 '+('0 '*9)+'\n')
  nonzero('negative Farkas weight',[WORK/'metric_audit15',q/'one.txt',q/'negative.classes',q/'residual.txt',q/'summary.txt'],expected='negative/missing multiplier')
 nonzero('disabled assertions in master',[sys.executable,'-O','-S','-B',ROOT/'verify_all.py','--help'],expected='requires assertions')
 nonzero('disabled assertions in force checker',[sys.executable,'-O','-S','-B',R/'exact_force15.py'],expected='without -O')
 out={'verified':True,'tests_passed':len(tests),'tests':tests,'seconds':time.time()-st}
 (R/'rejection15_verified.json').write_text(json.dumps(out,indent=2));print('REJECTION TESTS VERIFIED',len(tests),flush=True);return out
if __name__=='__main__':main()
