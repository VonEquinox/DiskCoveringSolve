#!/usr/bin/env python3
"""Fail-closed, fresh, standard-library accepting pipeline for the n=18 theorem."""
from pathlib import Path
import sys, os, json, hashlib, subprocess, time
R=Path(__file__).resolve().parent
REPORT=R/'reports';REPORT.mkdir(exist_ok=True)
MASTER=REPORT/'MASTER_VERIFIED.json'
MASTER.unlink(missing_ok=True)
if not __debug__ or os.environ.get('PYTHONOPTIMIZE','') not in ('','0'):
    raise RuntimeError('Do not disable assertions (-O, -OO, PYTHONOPTIMIZE)')

def inputs():
    out={}
    for p in sorted(R.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(R)
        if any(x in ('reports','__pycache__','discovery') for x in rel.parts):continue
        if p.name.endswith('_verified.json') or p.suffix in ('.pyc','.log'):continue
        out[str(rel)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out

def main():
    start=time.time();before=inputs();steps=[]
    required=['PROOF_zh.md','README.md','core/system18.py','core/root18.json',
              'enumeration/peel_enum.cpp','enumeration/peel_audit.cpp',
              'core/candidate18_tree.json.gz','enumeration/noncandidate18_trees.jsonl.gz']
    assert all(x in before for x in required),'missing proof input'
    sequence=['verify_interfaces18.py','verify_root18.py','verify_upper18.py',
              'anchor_isolation18.py','verify_enumeration18.py','verify_forest18.py','test_fail_closed18.py']
    for name in sequence:
        step=time.time()
        proc=subprocess.run([sys.executable,'-S','-B',str(R/'core'/name)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,cwd=R)
        (REPORT/(name+'.log')).write_text(proc.stdout,encoding='utf8')
        if proc.returncode:raise RuntimeError(f'{name} failed (exit {proc.returncode})\n{proc.stdout}')
        steps.append({'step':name,'seconds':time.time()-step,'log_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest()})
        print('PASSED',name,flush=True)
    def read(name):return json.loads((R/'core'/name).read_text())
    root=read('root18_verified.json');upper=read('upper18_verified.json')
    local=read('anchor_isolation_verified.json');interfaces=read('interfaces18_verified.json')
    enumeration=read('enumeration18_verified.json');forest=read('forest18_verified.json');tests=read('negative18_verified.json')
    assert all(x['verified'] is True for x in (root,upper,local,interfaces,enumeration,forest,tests))
    rh=before['core/root18.json']
    assert all(x['root_sha256']==rh for x in (root,upper,local,interfaces,forest))
    assert enumeration['residual_sha256']==before['enumeration/metric18_residuals.json']
    assert forest['local_certificate_sha256']==before['core/anchor_isolation_certificate.json']
    assert forest['candidate']['certificate_sha256']==before['core/candidate18_tree.json.gz']
    assert forest['noncandidate']['certificate_sha256']==before['enumeration/noncandidate18_trees.jsonl.gz']
    assert forest['noncandidate']['cases']+1==enumeration['surviving_orbits']
    assert interfaces['anchor_radius']==local['anchor_Frobenius_radius']
    assert inputs()==before,'accepting inputs changed during replay'
    out={'verified':True,'theorem':'r_18 = sqrt(t_star), R_18 = 1/sqrt(t_star)',
         'global_optimality_certificate_chain_closed':True,'kkt_dimension':root['dimension'],
         'decimal_intervals':root['decimal_intervals'],'root_sha256':rh,
         'covering_upper_verified':True,'upper_cover':{k:upper[k] for k in ('vertices','edges','faces','boundary_caps','exact_equality_center_faces')},
         'enumeration_method':enumeration['method'],'surviving_topology_orbits':enumeration['surviving_orbits'],
         'primary_partial_states':enumeration['primary_states'],'independent_partial_states':enumeration['independent_states'],
         'families':[{'B':x['B'],'I':x['I'],'surviving_orbits':x['primary']['orbits'],
                      'primary_states':x['primary']['calls'],'independent_states':x['independent']['visits'],
                      'rooted_survivors':x['primary']['rooted_survivors'],'manifest_sha256':x['manifest_sha256']} for x in enumeration['families']],
         'candidate_tree':forest['candidate'],'noncandidate_forest':forest['noncandidate'],
         'unresolved_leaves':0,'anchor_isolation':local,
         'small_unpruned_count_tests':len(enumeration['unpruned_count_regressions']),
         'negative_tests':tests,'input_sha256':before,'steps':steps,'seconds':time.time()-start,
         'trust_boundary':'Written convex-geometric, planar-topological, and exhaustive-peeling proofs; Python integer/Fraction arithmetic and a C++17 compiler. Every input certificate is replayed. No floating-point optimizer is an accepting dependency. Not Lean/Coq formalization, third-party review, or independent peer review.'}
    MASTER.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf8')
    print('R18 COMPLETE CERTIFICATE CHAIN VERIFIED',flush=True)
    print(json.dumps({k:out[k] for k in ('verified','theorem','decimal_intervals','surviving_topology_orbits','primary_partial_states','independent_partial_states','unresolved_leaves','seconds')},indent=2),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException:
        MASTER.unlink(missing_ok=True)
        raise
