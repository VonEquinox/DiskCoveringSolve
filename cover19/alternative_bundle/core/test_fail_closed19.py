"""Run corrupt-certificate tests against the actual accepting routines."""
from contextlib import contextmanager,redirect_stdout
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import gzip,io,json,shutil,subprocess,sys,tempfile,time
import geometry19 as S
import exact_force19 as G
import anchor_isolation19 as A
import verify_peeling19 as P
import verify_forest19 as V
import verify_interfaces19 as I
from exact_arcs import sincos_turn,SCALE
if not __debug__:raise RuntimeError('Run without -O')
ROOT=Path(__file__).resolve().parents[1]
@contextmanager
def patch(obj,key,value):
 old=getattr(obj,key);setattr(obj,key,value)
 try:yield
 finally:setattr(obj,key,old)

def main():
 st=time.time();results={}
 def rejected(name,fn):
  try:
   with redirect_stdout(io.StringIO()):fn()
  except (AssertionError,ValueError,RuntimeError,KeyError,IndexError,TypeError):results[name]='rejected'
  else:raise AssertionError('Invalid input ACCEPTED: '+name)
 rows=json.load(open(ROOT/'enumeration'/'survivors19.json'));cm=json.load(open(ROOT/'enumeration'/'candidate19_map.json'));candidate=cm['residual_index']
 with gzip.open(ROOT/'enumeration'/'noncandidate19_trees.jsonl.gz','rt') as f:
  next(f);case=next(c for c in map(json.loads,f) if len(c['nodes'])==1 and c['nodes'][0]['kind']=='DUAL')
 B=case['B'];fs=case['faces'];base=case['nodes'][0]
 def check(node):V.traverse(B,fs,[node],False)
 for field in ('force_num','support_num'):
  node=deepcopy(base);node[field][0][0]+=1
  rejected('changed_integer_'+field,lambda node=node:check(node))
 node=deepcopy(base);node['lambda_num'][0]=-1;rejected('negative_multiplier',lambda:check(node))
 node=deepcopy(base);node['lambda_num'][0]=True;rejected('boolean_multiplier',lambda:check(node))
 node=deepcopy(base);node['force_num'].pop();rejected('missing_rod_force',lambda:check(node))
 node=deepcopy(base);node['force_num'][0][0]=float(node['force_num'][0][0]);rejected('floating_force_field',lambda:check(node))
 node=deepcopy(base);node['force_num']=[[0,0] for _ in node['force_num']];node['support_num']=[[0,0] for _ in node['support_num']];node['lambda_num']=[0]*len(node['lambda_num']);rejected('zero_force_no_separation',lambda:check(node))
 for name,node in [('unresolved_leaf',{'kind':'UNRESOLVED'}),('false_empty_leaf',{'kind':'EMPTY'}),('noncandidate_local_leaf',{'kind':'LOCAL'})]:rejected(name,lambda node=node:check(node))
 rejected('unreachable_node',lambda:V.traverse(B,fs,[deepcopy(base),deepcopy(base)],False))
 with gzip.open(ROOT/'core'/'candidate19_tree.json.gz','rt') as f:tree=json.load(f)
 assert tree['nodes'][0]['kind']=='SPLIT'
 for name,mut in [('duplicate_child',lambda n:n['children'].__setitem__(1,n['children'][0])),('cycle',lambda n:n['children'].__setitem__(0,0)),('child_out_of_range',lambda n:n['children'].__setitem__(0,10**8)),('invalid_axis',lambda n:n.__setitem__('axis',S.B)),('boolean_child',lambda n:n['children'].__setitem__(0,True)),('split_outside_parent',lambda n:n.__setitem__('mid','-1'))]:
  ns=deepcopy(tree['nodes']);mut(ns[0]);rejected(name,lambda ns=ns:V.traverse(12,S.FACES,ns,True))
 rejected('root_box_falsely_local',lambda:V.traverse(12,S.FACES,[{'kind':'LOCAL'}],True))
 complete=set(range(len(rows)))-{candidate}
 rejected('missing_noncandidate',lambda:V.require_complete_cases(complete-{0},len(rows),candidate))
 rejected('candidate_in_noncandidate',lambda:V.require_complete_cases(complete|{candidate},len(rows),candidate))
 bad=deepcopy(rows[0]['faces']);bad.append(bad[0]);rejected('duplicate_triangle',lambda:P.bfs_code(B,bad))
 bad=deepcopy(rows[0]['faces']);bad[0][-1]=19;rejected('triangle_vertex_out_of_range',lambda:P.bfs_code(B,bad))
 with patch(S,'ANCHORS',[(F(2),F(0))]+S.ANCHORS[1:]):rejected('corrupted_exact_anchor',S.check_geometry)
 with patch(S,'T',F(1,14)):rejected('wrong_exact_radius',S.check_geometry)
 ww=S.WEIGHTS.copy();ww[0]+=F(1,10**12)
 with patch(S,'WEIGHTS',ww):rejected('corrupted_exact_stress',S.check_geometry)
 with patch(S,'EDGES',S.EDGES[:-1]):rejected('missing_stressed_edge',I.verify)
 with patch(A,'RADIUS',F(1,2)):rejected('unproved_local_radius',A.prepare)
 with patch(G,'RU',F(1,4)):rejected('threshold_below_target',P.check_constants)
 with patch(G,'C6',1):rejected('invalid_six_rod_cap',P.check_constants)
 with tempfile.TemporaryDirectory(prefix='r19_reject_') as td:
  temp=Path(td);shutil.copyfile(ROOT/'core'/'geometry19.py',temp/'geometry19.py')
  lc=json.load(open(ROOT/'core'/'anchor_isolation_certificate.json'))
  bad=deepcopy(lc);bad['Rnum'][0][0]=0;(temp/'anchor_isolation_certificate.json').write_text(json.dumps(bad))
  with patch(A,'R',temp):rejected('singular_congruence',A.verify)
  bad=deepcopy(lc);bad['geometry_sha256']='0'*64;(temp/'anchor_isolation_certificate.json').write_text(json.dumps(bad))
  with patch(A,'R',temp):rejected('certificate_wrong_geometry',A.verify)
  bad=deepcopy(cm);bad['map'][0],bad['map'][1]=bad['map'][1],bad['map'][0];(temp/'candidate19_map.json').write_text(json.dumps(bad))
  with patch(V,'ENUM',temp):rejected('wrong_candidate_graph_mapping',lambda:V.check_mapping(rows))
  cmds=[('assertions_disabled',[sys.executable,'-S','-O',str(ROOT/'core'/'verify_upper19.py')]),('unsupported_dimension',[str(ROOT/'_runtime'/'peel_bfs19'),'20','12',str(temp/'unsupported.json')]),('bad_boundary_size',[str(ROOT/'_runtime'/'peel_bfs19'),'19','2',str(temp/'bad.json')]),('limited_search_rejected',[str(ROOT/'_runtime'/'peel_dfs19'),'19','12',str(temp/'limit.json'),'1'])]
  for name,cmd in cmds:
   p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
   assert p.returncode!=0,(name,p.stdout);results[name]='rejected'
  assert not(temp/'limit.json').exists()
 co,si=sincos_turn(F(1,12));assert F(si[0],SCALE)<=F(1,2)<=F(si[1],SCALE)
 co,si=sincos_turn(F(1,6));assert F(co[0],SCALE)<=F(1,2)<=F(co[1],SCALE)
 co,si=sincos_turn(F(1,4));assert co[0]<=0<=co[1] and si[0]<=SCALE<=si[1]
 out={'verified':True,'tests':results,'rejection_tests':len(results),'analytic_trig_checks':3,'seconds':time.time()-st}
 (ROOT/'reports'/'fail_closed19_verified.json').write_text(json.dumps(out,indent=2));print('FAIL-CLOSED TESTS VERIFIED',len(results),'rejections and 3 exact trig identities',flush=True);return out
if __name__=='__main__':main()
