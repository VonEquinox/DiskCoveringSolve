#!/usr/bin/env python3
"""Fresh, fail-closed replay of every n=18 proof component.

Run: python3 -S -B verify_all.py
No numerical optimizer, third-party Python package, stored success flag, or
previously generated topology output is accepted in place of a fresh check.
"""
from pathlib import Path
from fractions import Fraction
import hashlib, json, os, shutil, subprocess, sys, time
if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
    raise RuntimeError('Run without -O/-OO/PYTHONOPTIMIZE')
ROOT=Path(__file__).resolve().parent
REPORTS=ROOT/'reports'; RUNTIME=ROOT/'_runtime'; CORE=ROOT/'core'
REPORTS.mkdir(exist_ok=True); RUNTIME.mkdir(exist_ok=True)
MASTER=REPORTS/'MASTER_VERIFIED.json'
MASTER.unlink(missing_ok=True)

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def input_paths():
    base=[ROOT/'verify_all.py',ROOT/'PROOF_zh.md',ROOT/'README.md']
    base+=list(CORE.glob('*.py'))
    base += [CORE/name for name in ['root18.json','anchor_angles18.json',
             'anchor_isolation_certificate.json','candidate18_tree.json.gz']]
    base+=list((ROOT/'enumeration').glob('*.cpp'))
    base += [ROOT/'enumeration'/name for name in ['survivors18.json',
             'candidate18_map.json','noncandidate18_trees.jsonl.gz']]
    assert all(p.is_file() for p in base)
    return sorted(set(base))

def replay():
    start=time.time();paths=input_paths();before={str(p.relative_to(ROOT)):digest(p) for p in paths}
    compiler=shutil.which('g++') or shutil.which('clang++')
    if compiler is None:raise RuntimeError('A C++17 compiler (g++ or clang++) is required')
    env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONPATH',None)
    steps=[]
    with (REPORTS/'FULL_REPLAY.log').open('w',encoding='utf8',buffering=1) as log:
        def emit(s):
            print(s,flush=True);log.write(s+'\n')
        def run(name,args,outputs=()):
            for p in outputs:p.unlink(missing_ok=True)
            t=time.time();emit('=== '+name)
            proc=subprocess.run(args,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            text=proc.stdout
            (REPORTS/(name+'.log')).write_text(text,encoding='utf8')
            log.write(text);log.flush()
            assert proc.returncode==0,(name,proc.returncode,text[-5000:])
            assert all(p.is_file() for p in outputs),(name,'missing output')
            item={'stage':name,'seconds':time.time()-t,'log_sha256':digest(REPORTS/(name+'.log'))}
            steps.append(item);emit('PASSED '+name)
        def py(name,where='core',outfile=None):
            args=[sys.executable,'-S','-B',str(ROOT/where/(name+'.py'))]
            outputs=[outfile] if outfile else []
            run(name,args,outputs)
            if outfile and outfile.parent==CORE:
                shutil.copyfile(outfile,REPORTS/outfile.name)
        py('verify_interfaces18',outfile=CORE/'interfaces18_verified.json')
        py('verify_root18',outfile=CORE/'root18_verified.json')
        py('verify_upper18',outfile=CORE/'upper18_verified.json')
        py('anchor_isolation18',outfile=CORE/'anchor_isolation_verified.json')
        for mode in ['dfs','bfs']:
            exe=RUNTIME/f'peel_{mode}18'
            run('compile_'+mode,[compiler,'-O3','-std=c++17',str(ROOT/'enumeration'/f'peel_{mode}18.cpp'),'-o',str(exe)],[exe])
        for b in range(11,18):
            for mode in ['dfs','bfs']:
                out=RUNTIME/f'{mode}_B{b}.json'
                run(f'{mode}_B{b}',[str(RUNTIME/f'peel_{mode}18'),'18',str(b),str(out)],[out])
        py('verify_peeling18',outfile=REPORTS/'peeling18_verified.json')
        py('verify_forest18',outfile=CORE/'forest18_verified.json')
        py('verify_recursion_tests18',outfile=REPORTS/'recursion_tests18_verified.json')
        py('test_fail_closed18',outfile=REPORTS/'fail_closed18_verified.json')
        def read(name):return json.loads((REPORTS/name).read_text())
        root=read('root18_verified.json');upper=read('upper18_verified.json')
        local=read('anchor_isolation_verified.json');part=read('peeling18_verified.json')
        forest=read('forest18_verified.json');negative=read('fail_closed18_verified.json')
        tests=read('recursion_tests18_verified.json');interfaces=read('interfaces18_verified.json')
        assert all(d.get('verified') is True for d in [root,upper,local,part,forest,negative,tests,interfaces])
        rh=digest(CORE/'root18.json')
        assert interfaces['root_sha256']==root['root_sha256']==upper['root_sha256']==local['root_sha256']==part['constants']['root_sha256']==forest['root_sha256']==rh
        assert forest['local_certificate_sha256']==local['certificate_sha256']==digest(CORE/'anchor_isolation_certificate.json')
        assert part['survivors_sha256']==digest(ROOT/'enumeration'/'survivors18.json')
        assert root['dimension']==174 and (upper['vertices'],upper['edges'],upper['faces'],upper['boundary_caps'])==(30,75,46,12)
        assert Fraction(local['anchor_Frobenius_radius'])==Fraction(11,100)
        assert Fraction(local['isolation_margin'])>0
        assert part['all_families']==list(range(11,18))
        assert part['surviving_orbits']==forest['noncandidate']['cases']+1
        ca=forest['candidate'];nc=forest['noncandidate']
        assert set(ca['counts'])=={'SPLIT','DUAL','LOCAL'} and set(nc['counts'])=={'SPLIT','DUAL'}
        assert sum(ca['counts'].values())==ca['nodes']==2*ca['counts']['SPLIT']+1
        assert sum(nc['counts'].values())==nc['nodes']==2*nc['counts']['SPLIT']+nc['cases']
        assert Fraction(ca['minimum_force_margin'])>0 and Fraction(ca['minimum_local_margin'])>0 and Fraction(nc['minimum_force_margin'])>0
        after={str(p.relative_to(ROOT)):digest(p) for p in paths}
        assert before==after,'An accepting input changed during replay'
        out={'verified':True,'theorem':'r_18 = sqrt(t_star), R_18 = 1/sqrt(t_star)',
             'global_optimality_certificate_chain_closed':True,'kkt_dimension':root['dimension'],
             'decimal_intervals':root['decimal_intervals'],'root_sha256':rh,
             'covering_upper_verified':True,'upper_complex':{'vertices':30,'edges':75,'triangles':46,'circular_caps':12},
             'exhaustive_partial_topology_search_closed':True,
             'topology_method':'Complete root-face peeling with exact partial-graph negative-cycle cuts; independent DFS/Bellman and BFS/Floyd re-execution and exact survivor-set equality.',
             'surviving_orbits':part['surviving_orbits'],'families':part['families'],
             'noncandidate_cases':nc['cases'],'candidate_cases':1,
             'candidate_tree':ca,'noncandidate_forest':nc,
             'unresolved_combinatorial_branches':0,'unresolved_leaves':0,
             'anchor_isolation_radius':local['anchor_Frobenius_radius'],
             'negative_tests':negative,'unpruned_recurrence_tests':tests,
             'input_sha256':before,'steps':steps,'seconds':time.time()-start,
             'trust_boundary':'Written geometric/topological and finite-recursion proofs in PROOF_zh.md; Python integer/Fraction arithmetic and standard C++17 compiler/runtime. No numerical optimizer or existing success report is trusted. Not Lean/Coq and not external peer review.'}
        tmp=REPORTS/'MASTER_VERIFIED.json.tmp';tmp.write_text(json.dumps(out,indent=2));tmp.replace(MASTER)
        emit('R18 COMPLETE CERTIFICATE CHAIN VERIFIED')
        emit(f"{nc['cases']} + 1 = {part['surviving_orbits']} surviving orbits after exhaustive partial cuts; unresolved leaves = 0")
        emit('All accepting input hashes unchanged during fresh replay.')
    return out

if __name__=='__main__':
    try:replay()
    except BaseException:
        MASTER.unlink(missing_ok=True)
        (REPORTS/'MASTER_VERIFIED.json.tmp').unlink(missing_ok=True)
        raise
