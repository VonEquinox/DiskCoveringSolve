"""Rejection tests against the actual accepting routines, standard library only."""
from contextlib import contextmanager, redirect_stdout
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import gzip, io, json, shutil, subprocess, sys, tempfile, time
import system18 as S
import exact_force18 as G
import verify_root18 as VR
import anchor_isolation18 as A
import verify_peeling18 as P
import verify_forest18 as V
from exact_arcs import sincos_turn,SCALE
if not __debug__:raise RuntimeError('Run without -O')
ROOT=Path(__file__).resolve().parents[1]

@contextmanager
def patch(obj,key,value):
    old=getattr(obj,key);setattr(obj,key,value)
    try:yield
    finally:setattr(obj,key,old)

def main():
    started=time.time();results={}
    def rejected(name,fn):
        try:
            with redirect_stdout(io.StringIO()):fn()
        except (AssertionError,ValueError,RuntimeError,KeyError,IndexError,TypeError):
            results[name]='rejected'
        else:
            raise AssertionError('Invalid certificate ACCEPTED: '+name)
    rows=json.loads((ROOT/'enumeration'/'survivors18.json').read_text())
    cmap=json.loads((ROOT/'enumeration'/'candidate18_map.json').read_text())
    cand=cmap['residual_index']
    with gzip.open(ROOT/'enumeration'/'noncandidate18_trees.jsonl.gz','rt') as f:
        json.loads(next(f));case=next(c for c in (json.loads(line) for line in f) if len(c['nodes'])==1 and c['nodes'][0]['kind']=='DUAL')
    assert len(case['nodes'])==1 and case['nodes'][0]['kind']=='DUAL'
    B=case['B'];faces=case['faces'];base=case['nodes'][0]
    def check_node(node):V.traverse(B,faces,[node],False)
    for field,position in [('force_num',(0,0)),('support_num',(0,0))]:
        node=deepcopy(base);node[field][position[0]][position[1]]+=1
        rejected('one_integer_'+field+'_unit_corrupted',lambda node=node:check_node(node))
    node=deepcopy(base);node['lambda_num'][0]=-1
    rejected('negative_halfplane_multiplier',lambda:check_node(node))
    node=deepcopy(base);node['lambda_num'][0]=True
    rejected('boolean_in_integer_multiplier',lambda:check_node(node))
    node=deepcopy(base);node['force_num'].pop()
    rejected('missing_edge_force',lambda:check_node(node))
    rejected('unresolved_leaf',lambda:check_node({'kind':'UNRESOLVED'}))
    rejected('false_empty_leaf',lambda:check_node({'kind':'EMPTY'}))
    rejected('local_leaf_in_noncandidate',lambda:check_node({'kind':'LOCAL'}))
    rejected('unreachable_extra_node',lambda:V.traverse(B,faces,[deepcopy(base),deepcopy(base)],False))
    with gzip.open(ROOT/'core'/'candidate18_tree.json.gz','rt') as f:tree=json.load(f)
    assert tree['nodes'][0]['kind']=='SPLIT'
    for name,mutate in [
        ('duplicate_split_child',lambda z:z['children'].__setitem__(1,z['children'][0])),
        ('cyclic_child',lambda z:z['children'].__setitem__(0,0)),
        ('out_of_range_child',lambda z:z['children'].__setitem__(0,10**9)),
        ('invalid_split_axis',lambda z:z.__setitem__('axis',S.B)),
        ('boolean_child_index',lambda z:z['children'].__setitem__(0,True)),
        ('cut_outside_parent',lambda z:z.__setitem__('mid','-1'))]:
        t=deepcopy(tree['nodes']);mutate(t[0])
        rejected(name,lambda t=t:V.traverse(S.B,S.FACES,t,True))
    rejected('whole_root_box_falsely_local',lambda:V.traverse(S.B,S.FACES,[{'kind':'LOCAL'}],True))
    complete=set(range(len(rows)))-{cand}
    rejected('missing_noncandidate_case',lambda:V.require_complete_cases(complete-{0},len(rows),cand))
    rejected('candidate_in_noncandidate_forest',lambda:V.require_complete_cases(complete|{cand},len(rows),cand))
    bad=deepcopy(rows[0]['faces']);bad.append(bad[0])
    rejected('duplicate_triangle',lambda:P.bfs_code(rows[0]['B'],bad))
    bad=deepcopy(rows[0]['faces']);bad[0][-1]=18
    rejected('out_of_range_triangle_vertex',lambda:P.bfs_code(rows[0]['B'],bad))
    with tempfile.TemporaryDirectory(prefix='r18-reject-') as td:
        temp=Path(td)
        root=json.loads((S.R/'root18.json').read_text())
        bad=deepcopy(root);bad['xnum'][S.TID]=str(int(bad['xnum'][S.TID])+10**40)
        (temp/'root18.json').write_text(json.dumps(bad))
        with patch(S,'R',temp),patch(VR,'R',temp):
            rejected('corrupted_algebraic_radius',VR.verify)
        bad=deepcopy(root);bad['Ynum'][0]=['0']*S.D
        (temp/'root18.json').write_text(json.dumps(bad))
        with patch(S,'R',temp),patch(VR,'R',temp):
            rejected('singular_root_preconditioner',VR.verify)
        shutil.copyfile(ROOT/'core'/'root18.json',temp/'root18.json')
        shutil.copyfile(ROOT/'core'/'system18.py',temp/'system18.py')
        lc=json.loads((ROOT/'core'/'anchor_isolation_certificate.json').read_text())
        bad=deepcopy(lc);bad['Rnum'][0][0]=0
        (temp/'anchor_isolation_certificate.json').write_text(json.dumps(bad))
        with patch(A,'ROOT',temp):rejected('singular_positive_matrix_congruence',A.main)
        bad=deepcopy(lc);bad['root_sha256']='0'*64
        (temp/'anchor_isolation_certificate.json').write_text(json.dumps(bad))
        with patch(A,'ROOT',temp):rejected('local_certificate_wrong_root',A.main)
        bad=deepcopy(cmap);bad['map'][0],bad['map'][1]=bad['map'][1],bad['map'][0]
        (temp/'candidate18_map.json').write_text(json.dumps(bad))
        with patch(V,'ENUM',temp):rejected('invalid_candidate_graph_mapping',lambda:V.check_mapping(rows))
        for name,cmd in [
            ('assertions_disabled',[sys.executable,'-S','-O',str(ROOT/'core'/'verify_root18.py')]),
            ('unsupported_combinatorial_size',[str(ROOT/'_runtime'/'peel_bfs18'),'19','11',str(temp/'invalid.json')]),
            ('invalid_boundary_size',[str(ROOT/'_runtime'/'peel_bfs18'),'18','2',str(temp/'invalid2.json')]),
            ('discovery_test_limit_not_complete',[str(ROOT/'_runtime'/'peel_dfs18'),'18','11',str(temp/'limited.json'),'1'])]:
            proc=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            assert proc.returncode!=0,(name,proc.stdout)
            results[name]='rejected'
        assert not (temp/'limited.json').exists()
    with patch(A,'RADIUS',F(1,2)):
        rejected('unsupported_local_radius',A.prepare)
    with patch(G,'RU',F(1,4)):
        rejected('threshold_below_candidate',P.check_constants)
    with patch(S,'EDGES',S.EDGES[:-1]):
        rejected('root_graph_edge_missing',VR.verify)
    # Analytic identities independently test the quadrant and pi conventions.
    co,si=sincos_turn(F(1,12));assert F(si[0],SCALE)<=F(1,2)<=F(si[1],SCALE)
    co,si=sincos_turn(F(1,6));assert F(co[0],SCALE)<=F(1,2)<=F(co[1],SCALE)
    co,si=sincos_turn(F(1,4));assert co[0]<=0<=co[1] and si[0]<=SCALE<=si[1]
    out={'verified':True,'tests':results,'rejection_tests':len(results),'analytic_trig_checks':3,'seconds':time.time()-started}
    (ROOT/'reports'/'fail_closed18_verified.json').write_text(json.dumps(out,indent=2))
    print('FAIL-CLOSED TESTS VERIFIED',len(results),'rejections and 3 exact trig identities',flush=True)
    return out
if __name__=='__main__':main()
