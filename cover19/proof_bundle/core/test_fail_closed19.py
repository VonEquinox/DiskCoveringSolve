"""Adversarial regression tests. These supplement, not replace, the proof."""
from pathlib import Path
from fractions import Fraction as F
import copy,json,gzip,time,tempfile,subprocess,sys,shutil
import exact_force19 as G
import geometry19 as S
import verify_forest19 as V
import verify_candidate19 as C
import verify_interfaces19 as I
from verify_enumeration19 import compile_tools,run
if not __debug__:raise RuntimeError('Assertions must be enabled')
R=Path(__file__).resolve().parent;E=R.parent/'enumeration'

def verify():
 (R/'negative19_verified.json').unlink(missing_ok=True);t=time.time();tests=[]
 def rejected(name,fn):
  try:fn()
  except (AssertionError,RuntimeError,ValueError,KeyError,IndexError,TypeError):tests.append(name);return
  raise RuntimeError('Corruption was accepted: '+name)
 def failed_process(name,args,**kwargs):
  p=subprocess.run(list(map(str,args)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,**kwargs)
  assert p.returncode!=0,(name,p.stdout);tests.append(name)
 rows=json.load(open(E/'metric19_residuals.json'));candidate=V.check_mapping(rows)
 with gzip.open(E/'noncandidate19_trees.jsonl.gz','rt') as f:head=json.loads(f.readline());case=json.loads(f.readline())
 assert len(case['nodes'])==1 and case['nodes'][0]['kind']=='DUAL'
 b=case['B'];faces=case['faces'];nodes=case['nodes'];ed=G.edges(b,faces);M=G.bounds(b,faces);lo,hi=G.tighten(M,*G.rootbox(M));cert=nodes[0]
 assert G.check(b,ed,lo,hi,cert)>0
 for name,key,action in [
 ('one integer force unit changed','force_num',lambda v:v[0].__setitem__(0,v[0][0]+1)),
 ('one integer support unit changed','support_num',lambda v:v[0].__setitem__(0,v[0][0]+1)),
 ('negative halfplane multiplier','lambda_num',lambda v:v.__setitem__(0,-1)),
 ('boolean multiplier','lambda_num',lambda v:v.__setitem__(0,True)),
 ('missing edge force','force_num',lambda v:v.pop())]:
  bad=copy.deepcopy(cert);action(bad[key]);rejected(name,lambda bad=bad:G.check(b,ed,lo,hi,bad))
 for name,kind in [('unresolved leaf','UNRESOLVED'),('false empty leaf','EMPTY'),('local leaf forbidden everywhere','LOCAL')]:
  rejected(name,lambda kind=kind:V.traverse(b,faces,[{'kind':kind}]))
 bad=copy.deepcopy(nodes)+[{'kind':'EMPTY'}];rejected('unreachable appended node',lambda:V.traverse(b,faces,bad))
 base={'kind':'SPLIT','axis':0,'mid':str((lo[0]+hi[0])/2),'children':[1,2]}
 changes=[('duplicate branch children','children',[1,1]),('cyclic branch','children',[0,1]),('out of range child','children',[20,1]),('boolean child index','children',[True,2]),('invalid split axis','axis',b),('boundary split point','mid',str(lo[0])),('outside split point','mid',str(lo[0]-1))]
 for name,key,val in changes:
  root=copy.deepcopy(base);root[key]=val;bad=[root,copy.deepcopy(cert),copy.deepcopy(cert)]
  rejected(name,lambda bad=bad:V.traverse(b,faces,bad))
 seen=set(range(len(rows)))-{candidate};badseen=seen-{next(iter(seen))}
 rejected('missing residual topology',lambda:V.complete_case_set(badseen,len(rows),candidate))
 rejected('candidate in noncandidate set',lambda:V.complete_case_set(seen|{candidate},len(rows),candidate))
 cm=json.load(open(E/'candidate19_map.json'));cm['map'][0],cm['map'][-1]=cm['map'][-1],cm['map'][0]
 rejected('wrong candidate graph map',lambda:V.check_mapping(rows,cm))
 rejected('nonpositive Hessian pivot',lambda:C.ldl([[F(-1)]]))
 with tempfile.TemporaryDirectory(prefix='cover19_negative_') as td:
  tmp=Path(td);oldcr,oldir=C.R,I.R;C.R=tmp;I.R=tmp
  try:
   for name,mod,key,val,fn in [
    ('wrong algebraic radius',S,'T',F(2,13),C.verify),
    ('missing positive rod',S,'EDGES',S.EDGES[:-1],C.verify),
    ('missing positive weight',S,'WEIGHTS',S.WEIGHTS[:-1],C.verify),
    ('incorrect anchor',S,'Q',[(F(0),F(1))]+S.Q[1:],C.verify),
    ('unsafe angular cap',G,'A',1,I.verify),
    ('threshold below candidate',G,'RU',F(1,4),I.verify)]:
    old=getattr(mod,key);setattr(mod,key,val)
    try:rejected(name,fn)
    finally:setattr(mod,key,old)
  finally:C.R=oldcr;I.R=oldir
  primary,audit=compile_tools(tmp)
  pfx=tmp/'small';report=tmp/'audit.json';run([primary,6,4,G.A,G.C,G.C6,pfx,0]);archive=pfx.with_suffix('.txt');content=archive.read_text()
  run([audit,6,4,G.A,G.C,G.C6,archive,report,0])
  for name,contents in [('duplicate topology representative',content+content.splitlines()[0]+'\n'),('out of range triangle mask','1 999999999999999999\n'),('truncated triangle record','1 7\n'),('noninteger triangle token','1 junk\n')]:
   bad=tmp/'bad.txt';bad.write_text(contents);report.write_text('{"verified":true}')
   failed_process(name,[audit,6,4,G.A,G.C,G.C6,bad,report,0]);assert not report.exists()
  pfx=tmp/'limited';pfx.with_suffix('.json').write_text('{"complete":true}');pfx.with_suffix('.txt').write_text('stale')
  failed_process('test-limit truncation',[primary,8,3,G.A,G.C,G.C6,pfx,0,1]);assert not pfx.with_suffix('.json').exists() and not pfx.with_suffix('.txt').exists()
  pfx.with_suffix('.json').write_text('{"complete":true}')
  failed_process('unsupported size',[primary,21,12,G.A,G.C,G.C6,pfx,1]);assert not pfx.with_suffix('.json').exists()
  failed_process('invalid boundary size',[primary,19,2,G.A,G.C,G.C6,pfx,1])
  failed_process('assertions disabled in checker',[sys.executable,'-S','-O','-c',f'import sys;sys.path.insert(0,{str(R)!r});import exact_force19'])
  failed_process('assertions disabled in master',[sys.executable,'-S','-O',R.parent/'verify_all.py'])
 out={'verified':True,'rejection_tests':len(tests),'tests':tests,'seconds':time.time()-t}
 (R/'negative19_verified.json').write_text(json.dumps(out,indent=2));print('FAIL-CLOSED TESTS VERIFIED',len(tests),flush=True);return out
if __name__=='__main__':verify()
