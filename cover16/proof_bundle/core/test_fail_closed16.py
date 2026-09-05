"""Adversarial regression tests; these supplement, not replace, exact replay."""
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import copy,gzip,json,shutil,struct,subprocess,sys,tempfile,time
import exact_force16 as G
import system16 as S
from verify_forest16 import traverse,check_mapping
from verify_partition16 import verify_lp,labelled
R=Path(__file__).resolve().parent;P=R.parent;E=P/'enumeration'

def rejected(fn):
 try:fn()
 except (AssertionError,ValueError,KeyError,IndexError,RuntimeError,StopIteration,TypeError):return True
 return False

def main():
 st=time.time();tests={}
 with gzip.open(E/'noncandidate16_trees.jsonl.gz','rt',encoding='utf8') as f:
  header=json.loads(f.readline());case=None
  for line in f:
   candidate=json.loads(line)
   if len(candidate['nodes'])==1 and candidate['nodes'][0]['kind']=='DUAL':case=candidate;break
 assert case is not None
 B=case['B'];faces=case['faces'];traverse(B,faces,case['nodes'],False)
 def reject_tree(name,nodes,allow=False,bb=B,ff=faces):
  assert rejected(lambda:traverse(bb,ff,nodes,allow)),name;tests[name]='rejected'
 bad=copy.deepcopy(case['nodes']);bad[0]['force_num'][0][0]+=1
 reject_tree('one_integer_force_unit_corrupted',bad)
 bad=copy.deepcopy(case['nodes']);bad[0]['lambda_num'][0]=-1
 reject_tree('negative_halfplane_multiplier',bad)
 bad=copy.deepcopy(case['nodes']);bad[0]['force_num'].pop()
 reject_tree('missing_rod_force',bad)
 reject_tree('unresolved_leaf',[{'kind':'UNRESOLVED'}])
 reject_tree('unreachable_extra_node',copy.deepcopy(case['nodes'])+[{'kind':'EMPTY'}])
 M=G.bounds(B,faces);lo,hi=G.rootbox(M);lo,hi=G.tighten(M,lo,hi);axis=next(i for i,(a,b) in enumerate(zip(lo,hi)) if a<b);mid=(lo[axis]+hi[axis])/2
 reject_tree('duplicate_split_child',[{'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[1,1]},{'kind':'EMPTY'}])
 reject_tree('cyclic_split_edge',[{'kind':'SPLIT','axis':axis,'mid':str(mid),'children':[0,1]},{'kind':'EMPTY'}])
 reject_tree('split_outside_parent',[{'kind':'SPLIT','axis':axis,'mid':str(lo[axis]-1),'children':[1,2]},{'kind':'EMPTY'},{'kind':'EMPTY'}])
 reject_tree('local_leaf_in_noncandidate_case',[{'kind':'LOCAL'}])
 reject_tree('whole_candidate_domain_falsely_local',[{'kind':'LOCAL'}],True,S.B,S.FACES)
 metric=json.load(open(E/'metric16_certificate.json'))
 fam=next(k for k,v in metric.items() if v['certified_signatures'])
 cb=int(fam.split('_')[0][1:]);ky,vals=next(iter(metric[fam]['certified_signatures'].items()));key=tuple(map(int,ky.split(','))) if ky else ()
 verify_lp(cb,key,vals);badvals=list(vals);badvals[0]=-1
 assert rejected(lambda:verify_lp(cb,key,badvals));tests['negative_Farkas_weight']='rejected'
 assert rejected(lambda:verify_lp(cb,key,[0]*len(vals)));tests['zero_Farkas_cover']='rejected'
 with tempfile.TemporaryDirectory(prefix='r16-reject-') as td:
  T=Path(td)
  for n in ['system16.py','verify_root16.py','root16.json','anchor_isolation16.py','anchor_isolation_certificate.json']:shutil.copyfile(R/n,T/n)
  c=json.load(open(T/'root16.json'));c['xnum'][0]=str(int(c['xnum'][0])+int(c['Qx'])//1000);(T/'root16.json').write_text(json.dumps(c))
  rr=subprocess.run([sys.executable,'-S','-B',str(T/'verify_root16.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'AssertionError' in rr.stdout;tests['corrupted_algebraic_root']='rejected'
  shutil.copyfile(R/'root16.json',T/'root16.json');c=json.load(open(T/'anchor_isolation_certificate.json'));c['Rnum'][0][0]=0;(T/'anchor_isolation_certificate.json').write_text(json.dumps(c))
  rr=subprocess.run([sys.executable,'-S','-B',str(T/'anchor_isolation16.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'AssertionError' in rr.stdout;tests['singular_PD_congruence']='rejected'
  with (P/'_runtime/B15_I1.masks').open('rb') as f:hh=list(struct.unpack('<5I',f.read(20)));rec=f.read(hh[3]*2)
  hh[4]=2;(T/'duplicate.masks').write_bytes(struct.pack('<5I',*hh)+rec+rec)
  (T/'empty.farkas').write_text('9974 21143 100000 1000000000000 0\n')
  rr=subprocess.run([str(P/'_runtime/audit16'),str(T/'duplicate.masks'),str(T/'empty.farkas'),str(labelled(15,1)),str(T/'duplicate_audit')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'duplicate orbit' in rr.stdout;tests['duplicate_orbit_representative']='rejected'
  hh[4]=1;(T/'invalid.masks').write_bytes(struct.pack('<5I',*hh)+b'\x00\x00'+rec[2:])
  rr=subprocess.run([str(P/'_runtime/audit16'),str(T/'invalid.masks'),str(T/'empty.farkas'),str(labelled(15,1)),str(T/'invalid_audit')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'bad triangle' in rr.stdout;tests['invalid_triangle_mask']='rejected'
  C=T/'missing_case'/'core';EN=T/'missing_case'/'enumeration';C.mkdir(parents=True);EN.mkdir()
  for n in ['system16.py','root16.json','interval16.py','exact_arcs.py','exact_force16.py','candidate_global16.py','verify_forest16.py','candidate16_tree.json.gz','anchor_isolation_certificate.json']:shutil.copyfile(R/n,C/n)
  for n in ['metric16_residuals.json','candidate16_map.json']:shutil.copyfile(E/n,EN/n)
  with gzip.open(EN/'noncandidate16_trees.jsonl.gz','wt',encoding='utf8') as f:f.write(json.dumps(header)+'\n'+json.dumps(case)+'\n')
  rr=subprocess.run([sys.executable,'-S','-B',str(C/'verify_forest16.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'AssertionError' in rr.stdout;tests['missing_noncandidate_cases']='rejected'
  # Candidate graph binding is tested without replaying the forest again.
  cm=json.load(open(EN/'candidate16_map.json'));cm['map'][0],cm['map'][1]=cm['map'][1],cm['map'][0];(EN/'candidate16_map.json').write_text(json.dumps(cm))
  code="import json;from pathlib import Path;import verify_forest16 as V;V.check_mapping(json.load(open(V.ENUM/'metric16_residuals.json')))"
  rr=subprocess.run([sys.executable,'-S','-B','-c',code],cwd=C,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert rr.returncode!=0 and b'AssertionError' in rr.stdout;tests['incorrect_candidate_graph_permutation']='rejected'
 rr=subprocess.run([sys.executable,'-S','-B','-O',str(P/'verify_all.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 assert rr.returncode!=0 and b'Run without' in rr.stdout;tests['assertions_disabled_in_master']='rejected'
 rr=subprocess.run([sys.executable,'-S','-B',str(R/'exact_force16.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 assert rr.returncode==0,'force checker must import normally before testing -O'
 rr=subprocess.run([sys.executable,'-S','-B','-O',str(R/'exact_force16.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 assert rr.returncode!=0 and b'without -O' in rr.stdout;tests['assertions_disabled_in_force_checker']='rejected'
 out={'verified':True,'tests_passed':len(tests),'tests':tests,'seconds':time.time()-st};(R/'fail_closed16_verified.json').write_text(json.dumps(out,indent=2));print('NEGATIVE TESTS PASSED',len(tests),tests,flush=True);return out
if __name__=='__main__':main()
