#!/usr/bin/env python3
"""Fresh replay of every part of the n=19 computer-assisted proof.

Run: python3 -S -B verify_all.py
Only Python's standard library and a C++17 compiler are needed. Numerical
proposal programs and pre-existing success reports are never proof inputs.
"""
from pathlib import Path
from fractions import Fraction
import hashlib, json, os, shutil, subprocess, sys, time
ROOT=Path(__file__).resolve().parent
REPORTS=ROOT/'reports'; RUNTIME=ROOT/'_runtime'; CORE=ROOT/'core'
REPORTS.mkdir(exist_ok=True); RUNTIME.mkdir(exist_ok=True)
MASTER=REPORTS/'MASTER_VERIFIED.json'
MASTER.unlink(missing_ok=True)
(REPORTS/'MASTER_VERIFIED.json.tmp').unlink(missing_ok=True)
if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
    raise RuntimeError('Run without -O/-OO/PYTHONOPTIMIZE')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def input_paths():
    paths=[ROOT/'verify_all.py',ROOT/'PROOF_zh.md',ROOT/'README.md']
    paths+=list(CORE.glob('*.py'))
    paths += [CORE/name for name in ['anchor_isolation_certificate.json','candidate19_tree.json.gz']]
    paths += [ROOT/'enumeration'/name for name in [
        'peel_dfs19.cpp','peel_bfs19.cpp','survivors19.json',
        'candidate19_map.json','noncandidate19_trees.jsonl.gz']]
    assert all(p.is_file() for p in paths), 'Missing accepting input'
    return sorted(set(paths))

def replay():
    started=time.time();paths=input_paths()
    before={str(p.relative_to(ROOT)):digest(p) for p in paths}
    compiler=shutil.which('g++') or shutil.which('clang++')
    if compiler is None: raise RuntimeError('A C++17 compiler (g++ or clang++) is required')
    env=dict(os.environ); env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONPATH',None)
    steps=[]
    with (REPORTS/'FULL_REPLAY.log').open('w',encoding='utf8',buffering=1) as log:
        def emit(s):
            print(s,flush=True);log.write(s+'\n')
        def run(name,args,outputs=()):
            for p in outputs:p.unlink(missing_ok=True)
            t=time.time();emit('=== '+name)
            proc=subprocess.run(args,cwd=ROOT,env=env,stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,text=True)
            text=proc.stdout
            lp=REPORTS/(name+'.log');lp.write_text(text,encoding='utf8')
            log.write(text);log.flush()
            assert proc.returncode==0,(name,proc.returncode,text[-5000:])
            assert all(p.is_file() for p in outputs),(name,'missing output')
            steps.append({'stage':name,'seconds':time.time()-t,'log_sha256':digest(lp)})
            emit('PASSED '+name)
        def py(name,outfile):
            run(name,[sys.executable,'-S','-B',str(CORE/(name+'.py'))],[outfile])
            if outfile.parent==CORE:shutil.copyfile(outfile,REPORTS/outfile.name)
        py('verify_interfaces19',CORE/'interfaces19_verified.json')
        py('verify_upper19',CORE/'upper19_verified.json')
        py('anchor_isolation19',CORE/'anchor_isolation_verified.json')
        for mode in ('dfs','bfs'):
            exe=RUNTIME/f'peel_{mode}19'
            run('compile_'+mode,[compiler,'-O3','-std=c++17',
                str(ROOT/'enumeration'/f'peel_{mode}19.cpp'),'-o',str(exe)],[exe])
        for b in range(12,19):
            for mode in ('dfs','bfs'):
                out=RUNTIME/f'{mode}_B{b}.json'
                run(f'{mode}_B{b}',[str(RUNTIME/f'peel_{mode}19'),'19',str(b),str(out)],[out])
        py('verify_peeling19',REPORTS/'peeling19_verified.json')
        py('verify_forest19',CORE/'forest19_verified.json')
        py('verify_recursion_tests19',REPORTS/'recursion_tests19_verified.json')
        py('test_fail_closed19',REPORTS/'fail_closed19_verified.json')
        def read(name):return json.loads((REPORTS/name).read_text())
        upper=read('upper19_verified.json');local=read('anchor_isolation_verified.json')
        part=read('peeling19_verified.json');forest=read('forest19_verified.json')
        negative=read('fail_closed19_verified.json');tests=read('recursion_tests19_verified.json')
        interfaces=read('interfaces19_verified.json')
        assert all(d.get('verified') is True for d in [upper,local,part,forest,negative,tests,interfaces])
        gh=digest(CORE/'geometry19.py')
        assert (interfaces['geometry_sha256']==upper['geometry_sha256']==local['geometry_sha256']
                ==part['constants']['geometry_sha256']==forest['geometry_sha256']==gh)
        assert (forest['local_certificate_sha256']==local['certificate_sha256']
                ==digest(CORE/'anchor_isolation_certificate.json'))
        assert part['survivors_sha256']==digest(ROOT/'enumeration'/'survivors19.json')
        assert Fraction(upper['t'])==Fraction(interfaces['radius_squared'])==Fraction(1,13)
        assert upper['exact_radius']=='1/sqrt(13)' and upper['exact_large_radius']=='sqrt(13)'
        assert (upper['vertices'],upper['edges'],upper['triangles'],upper['caps'])==(31,78,48,12)
        assert upper['degenerate_ring_faces']==[0,4,8,12,16,20]
        assert upper['center_equilateral_faces']==24
        assert (local['dimension'],local['rods'])==(70,42)
        assert interfaces['stressed_geometric_dimension']==70 and interfaces['stressed_rods']==42
        assert interfaces['unused_center_index']==12 and interfaces['full_auxiliary_nodes']==55
        assert Fraction(local['anchor_radius'])==Fraction(1,5)
        assert Fraction(local['isolation_margin'])==Fraction(17,76050)>0
        assert Fraction(local['congruence_margin'])>0
        assert part['all_families']==list(range(12,19)) and part['unresolved_combinatorial_branches']==0
        assert part['surviving_orbits']==forest['noncandidate']['cases']+1==24127
        assert part['rooted_survivors']==621582
        assert sum(d['surviving_orbits'] for d in part['families'].values())==24127
        ca=forest['candidate'];nc=forest['noncandidate']
        assert set(ca['counts'])=={'SPLIT','DUAL','LOCAL'} and set(nc['counts'])=={'SPLIT','DUAL'}
        assert sum(ca['counts'].values())==ca['nodes']==2*ca['counts']['SPLIT']+1
        assert sum(nc['counts'].values())==nc['nodes']==2*nc['counts']['SPLIT']+nc['cases']
        assert ca['nodes']==527 and ca['counts']=={'SPLIT':263,'DUAL':236,'LOCAL':28}
        assert nc['nodes']==31708 and nc['counts']=={'SPLIT':3791,'DUAL':27917}
        assert Fraction(ca['minimum_force_margin'])>0 and Fraction(ca['minimum_local_margin'])>0
        assert Fraction(nc['minimum_force_margin'])>0
        assert ca['certificate_sha256']==digest(CORE/'candidate19_tree.json.gz')
        assert nc['certificate_sha256']==digest(ROOT/'enumeration'/'noncandidate19_trees.jsonl.gz')
        assert negative['rejection_tests']==36 and tests['cases']==29
        after={str(p.relative_to(ROOT)):digest(p) for p in paths}
        assert before==after,'An accepting input changed during replay'
        out={'verified':True,'theorem':'r_19 = 1/sqrt(13), R_19 = sqrt(13)',
             'global_optimality_certificate_chain_closed':True,'exact_radius_squared':'1/13',
             'exact_candidate_field':'Q(sqrt(3))','root_isolation_required':False,
             'decimal_intervals':upper['decimal_intervals'],'geometry_sha256':gh,
             'covering_upper_verified':True,
             'upper_complex':{'vertices':31,'edges':78,'triangles':48,'circular_caps':12,
                              'exactly_degenerate_triangles':6},
             'positive_stress_graph':{'nodes':36,'rods':42,'geometric_dimension':70,
                                      'omitted_original_center':12},
             'exhaustive_partial_topology_search_closed':True,
             'topology_method':'Complete root-face peeling with exact partial-graph negative-cycle cuts; fresh DFS/Bellman and BFS/Floyd re-execution, independent normalization and exact survivor-set equality.',
             'surviving_orbits':part['surviving_orbits'],'rooted_survivors':part['rooted_survivors'],
             'families':part['families'],'noncandidate_cases':nc['cases'],'candidate_cases':1,
             'candidate_tree':ca,'noncandidate_forest':nc,
             'unresolved_combinatorial_branches':0,'unresolved_leaves':0,
             'anchor_isolation_radius':local['anchor_radius'],'anchor_isolation':local,
             'negative_tests':negative,'unpruned_recurrence_tests':tests,
             'input_sha256':before,'steps':steps,'seconds':time.time()-started,
             'trust_boundary':'Written geometric/topological and exhaustive finite-recursion proofs in PROOF_zh.md; Python integer/Fraction arithmetic and a standard C++17 compiler/runtime. No optimizer, floating eigenvalue, algebraic-root claim, or old success report is trusted. Not Lean/Coq and not external peer review.'}
        tmp=REPORTS/'MASTER_VERIFIED.json.tmp';tmp.write_text(json.dumps(out,indent=2));tmp.replace(MASTER)
        emit('R19 COMPLETE CERTIFICATE CHAIN VERIFIED')
        emit('r_19 = 1/sqrt(13); R_19 = sqrt(13).')
        emit(f"{nc['cases']} + 1 = {part['surviving_orbits']} surviving orbits after exhaustive partial cuts; unresolved leaves = 0")
        emit('All accepting input hashes unchanged during fresh replay.')
    return out

if __name__=='__main__':
    try:replay()
    except BaseException:
        MASTER.unlink(missing_ok=True)
        (REPORTS/'MASTER_VERIFIED.json.tmp').unlink(missing_ok=True)
        raise
