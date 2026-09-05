"""Deliberate corruption tests of the actual accepting routines.
These are regression tests, not replacements for the mathematical reduction.
"""
from pathlib import Path
from fractions import Fraction as F
import copy,gzip,json,subprocess,tempfile,struct,sys,time,shutil,random
from cases20 import Store,ROOT,ENUM
import domain20 as D
import compact_force20 as C
from verify_tree20 import verify as tree_verify
from verify_forest20 import consume_key,finish_keys
if not __debug__:raise RuntimeError('Run without -O')

def verify():
 st=time.time();store=Store();first=json.loads(next(iter(gzip.open(ENUM/'noncandidate20_trees.jsonl.gz','rt'))))
 B,j=first['B'],first['idx'];faces=store.faces(B,j);mat=D.bounds(B,faces);lo,hi=D.rootbox(mat);lo,hi=D.tighten(mat,lo,hi);ed=D.edges(B,faces)
 assert len(first['nodes'])==1 and first['nodes'][0]['kind']=='DUAL'
 cert=first['nodes'][0]['cert'];assert C.check(B,ed,lo,hi,cert)>0
 tests={}
 def reject(name,fn):
  try:fn()
  except (AssertionError,ValueError,RuntimeError,IndexError,KeyError,TypeError,OverflowError):tests[name]='rejected';return
  raise AssertionError(('corruption accepted',name))
 def mutated(key,value):
  c=copy.deepcopy(cert);c[key]=value;return lambda:C.check(B,ed,lo,hi,c)
 for value,label in [(0,'zero'),(-1,'negative'),(True,'boolean'),('1000000','string'),(10**19,'oversize')]:reject('force_denominator_'+label,mutated('q',value))
 a=cert['l'].copy();a[0]=-1;reject('negative_support_multiplier',mutated('l',a))
 a=cert['l'].copy();a[0]=1.5;reject('floating_support_multiplier',mutated('l',a))
 reject('missing_support_multiplier',mutated('l',cert['l'][:-1]))
 reject('missing_anchor_force_coordinate',mutated('b',cert['b'][:-1]))
 a=cert['b'].copy();a[0]=True;reject('boolean_anchor_force',mutated('b',a))
 reject('missing_cycle_force_coordinate',mutated('f',cert['f'][:-1]))
 a=cert['f'].copy();a[0]=1.5;reject('floating_cycle_force',mutated('f',a))
 a=copy.deepcopy(cert);a['unused']=0;reject('undeclared_force_field',lambda:C.check(B,ed,lo,hi,a))
 a=copy.deepcopy(cert);del a['q'];reject('missing_force_field',lambda:C.check(B,ed,lo,hi,a))
 zero={k:([0]*len(v) if isinstance(v,list) else v) for k,v in cert.items()}
 reject('zero_force_nonseparation',lambda:C.check(B,ed,lo,hi,zero))
 def tv(nodes):return lambda:tree_verify(B,faces,nodes,False)
 reject('unresolved_leaf',tv([{'kind':'UNRESOLVED'}]))
 reject('local_leaf_in_noncandidate',tv([{'kind':'LOCAL'}]))
 reject('false_empty_leaf',tv([{'kind':'EMPTY'}]))
 reject('unreachable_extra_node',tv(copy.deepcopy(first['nodes'])+[{'kind':'EMPTY'}]))
 reject('undeclared_node_field',tv([dict(first['nodes'][0],extra=1)]))
 split={'kind':'SPLIT','axis':0,'mid':str((lo[0]+hi[0])/2),'children':[1,2]}
 for name,key,value in [('duplicate_children','children',[1,1]),('missing_child','children',[1,9]),('cycle_to_parent','children',[0,1]),('boolean_axis','axis',True),('negative_axis','axis',-1),('invalid_axis','axis',B),('endpoint_split','mid',str(lo[0])),('outside_split','mid',str(hi[0]+1))]:
  a=copy.deepcopy(split);a[key]=value;reject(name,tv([a,{'kind':'EMPTY'},{'kind':'EMPTY'}]))
 reject('missing_child_field',tv([{'kind':'SPLIT','axis':0,'mid':str((lo[0]+hi[0])/2)}]))
 reject('duplicate_case',lambda:(consume_key(iter([(12,0)]),{'B':12,'idx':1})))
 reject('unexpected_extra_case',lambda:consume_key(iter([]),{'B':12,'idx':0}))
 reject('missing_final_case',lambda:finish_keys(iter([(12,0)])))
 reject('reordered_case',lambda:consume_key(iter([(13,4),(12,1)]),{'B':12,'idx':1}))
 # The integer-lattice tightening must equal direct Fraction propagation.
 rng=random.Random(2020)
 for _ in range(100):
  ll=[];hh=[]
  for l,h in zip(lo,hi):
   u=F(rng.randrange(129),128);v=F(rng.randrange(129),128);ll.append(l+(h-l)*min(u,v));hh.append(l+(h-l)*max(u,v))
  l=[F(0)]+ll;h=[F(0)]+hh
  lower=[max(l[i]-F(mat[k][i],D.D) for i in range(B)) for k in range(B)]
  upper=[min(h[i]+F(mat[i][k],D.D) for i in range(B)) for k in range(B)]
  ref=None if any(a>b for a,b in zip(lower,upper)) else (lower[1:],upper[1:])
  assert D.tighten(mat,ll,hh)==ref
 with tempfile.TemporaryDirectory(prefix='negative20_',dir=ROOT/'_runtime') as td:
  td=Path(td)
  def command_rejected(name,args):
   proc=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
   assert proc.returncode!=0,(name,proc.stdout,proc.stderr)
   tests[name]='rejected'
  # Compile-time and runtime guards are exercised through actual subprocesses.
  command_rejected('assertions_disabled',[sys.executable,'-O',str(ROOT/'verify_all.py')])
  guard=td/'guard';guard.mkdir();(guard/'reports').mkdir();shutil.copy(ROOT/'verify_all.py',guard/'verify_all.py')
  stale=guard/'reports/MASTER_VERIFIED.json';stale.write_text('{"verified":true}')
  command_rejected('disabled_assertions_remove_stale_master',[sys.executable,'-O',str(guard/'verify_all.py')])
  assert not stale.exists(),'stale master survived the rejected invocation'

  for n,b in [(21,12),(20,2),(20,21)]:
   command_rejected(f'invalid_recursion_dimensions_{n}_{b}',[str(ROOT/'_runtime/peel_bfs20'),str(n),str(b),str(td/'bad.json')])
  mask=list(store.masks(B,j));a=td/'a.txt';b=td/'b.txt';b.write_text('0\n')
  a.write_text('2\n'+' '.join(map(str,mask))+'\n'+' '.join(map(str,mask))+'\n')
  command_rejected('duplicate_orbit_representative',[str(ROOT/'_runtime/compare_maps20'),'20',str(B),str(a),str(b)])
  mask[0]=3;a.write_text('1\n'+' '.join(map(str,mask))+'\n')
  command_rejected('invalid_triangle_mask',[str(ROOT/'_runtime/compare_maps20'),'20',str(B),str(a),str(b)])
  copycore=td/'core';copycore.mkdir()
  for fn in ['system20.py','root20.json','verify_root20.py','anchor_isolation20.py','anchor_isolation_certificate.json']:shutil.copy(ROOT/'core'/fn,copycore/fn)
  root=json.load(open(copycore/'root20.json'));root['xnum'][0]=str(int(root['xnum'][0])+int(root['Qx']));(copycore/'root20.json').write_text(json.dumps(root))
  command_rejected('corrupted_algebraic_root',[sys.executable,'-S','-B',str(copycore/'verify_root20.py')])
  shutil.copy(ROOT/'core/root20.json',copycore/'root20.json')
  ac=json.load(open(copycore/'anchor_isolation_certificate.json'));ac['Rnum'][0][0]=0;(copycore/'anchor_isolation_certificate.json').write_text(json.dumps(ac))
  command_rejected('singular_congruence_matrix',[sys.executable,'-S','-B',str(copycore/'anchor_isolation20.py')])
  # A malformed packed file must fail before any topology can be used.
  import cases20 as CS
  old=CS.ENUM;CS.ENUM=td
  try:
   gzip.open(td/'survivors20_B12.bin.gz','wb').write(struct.pack('<4sIII',b'FAIL',20,12,0))
   reject('invalid_packed_manifest_header',lambda:CS.Store())
  finally:CS.ENUM=old
 out={'verified':True,'tests':tests,'negative_test_count':len(tests),'exact_integer_lattice_crosschecks':100,'seconds':time.time()-st}
 (ROOT/'reports/negative20_verified.json').write_text(json.dumps(out,indent=2));print('CORRUPTION TESTS REJECTED',len(tests),'LATTICE CROSSCHECKS',100,flush=True);return out
if __name__=='__main__':verify()
