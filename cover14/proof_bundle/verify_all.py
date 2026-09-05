#!/usr/bin/env python3
"""Fresh, fail-closed replay of the entire r14 proof certificate chain."""
if not __debug__:
    raise RuntimeError('Exact verification must not run with -O/-OO/PYTHONOPTIMIZE')
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import gzip,hashlib,json,os,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parent;CORE=ROOT/'core';ENUM=ROOT/'enumeration';REPORT=ROOT/'reports';RUNTIME=ROOT/'_runtime'
MASTER=REPORT/'MASTER_VERIFIED.json'
INPUTS=[
 'core/system14.py','core/root14.json','core/verify_root14.py','core/symmetry14.py','core/verify_upper14.py','core/anchor_angles14.json',
 'core/anchor_isolation14.py','core/anchor_isolation_certificate.json','core/exact_arcs.py','core/exact_force14.py','core/candidate_global14.py',
 'core/candidate14_tree.json.gz','core/verify_forest14.py','core/verify_partition14.py','core/test_fail_closed14.py',
 'enumeration/audit14.cpp','enumeration/metric14_certificate.json','enumeration/metric14_residuals.json','enumeration/candidate14_map.json',
 'enumeration/noncandidate14_trees.json.gz',*['enumeration/B%d_I%d.txt.gz'%(B,14-B) for B in range(10,14)],'verify_all.py']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run_python(name,outputs):
    for n in outputs:(CORE/n).unlink(missing_ok=True)
    log=REPORT/(name+'.log');st=time.time()
    with log.open('wb') as f:subprocess.run([sys.executable,'-S','-B',str(CORE/name)],cwd=CORE,stdout=f,stderr=subprocess.STDOUT,check=True)
    results={}
    for n in outputs:
        p=CORE/n;d=json.loads(p.read_text());assert d['verified'] is True;results[n]=d;shutil.move(p,REPORT/n)
    print('PASSED',name,flush=True)
    return results,{'step':name,'seconds':time.time()-st,'log_sha256':sha(log)}

def audit_family(B):
    p=RUNTIME/f'B{B}_I{14-B}.txt';out=REPORT/f'audit{B}.json';out.unlink(missing_ok=True);log=REPORT/f'audit{B}.log';st=time.time()
    with gzip.open(ENUM/(p.name+'.gz'),'rb') as f,p.open('wb') as g:shutil.copyfileobj(f,g)
    with log.open('wb') as f:subprocess.run([str(RUNTIME/'audit14'),str(p),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True)
    d=json.loads(out.read_text());assert d['verified'] is True;p.unlink()
    print('PASSED independent orbit audit',B,14-B,flush=True)
    return d,{'step':f'orbit_audit_{B}','seconds':time.time()-st,'log_sha256':sha(log)}

def main():
    REPORT.mkdir(exist_ok=True);RUNTIME.mkdir(exist_ok=True);MASTER.unlink(missing_ok=True);st=time.time();steps=[]
    initial={n:sha(ROOT/n) for n in INPUTS}
    results,s=run_python('verify_root14.py',['root14_verified.json']);root=results['root14_verified.json'];steps.append(s)
    results,s=run_python('verify_upper14.py',['upper14_verified.json','symmetry14_verified.json']);upper=results['upper14_verified.json'];sym=results['symmetry14_verified.json'];steps.append(s)
    results,s=run_python('anchor_isolation14.py',['anchor_isolation_verified.json']);local=results['anchor_isolation_verified.json'];steps.append(s)
    compiler=shutil.which('g++') or shutil.which('clang++')
    if compiler is None:raise RuntimeError('A C++17 compiler and Boost.Multiprecision headers are required')
    with (REPORT/'compile.log').open('wb') as f:subprocess.run([compiler,'-O3','-std=c++17',str(ENUM/'audit14.cpp'),'-o',str(RUNTIME/'audit14')],stdout=f,stderr=subprocess.STDOUT,check=True)
    with ThreadPoolExecutor(max_workers=4) as ex:
        audits=list(ex.map(audit_family,range(10,14)))
    for d,s in audits:steps.append(s)
    results,s=run_python('verify_partition14.py',['partition14_verified.json']);partition=results['partition14_verified.json'];steps.append(s)
    results,s=run_python('verify_forest14.py',['forest14_verified.json']);forest=results['forest14_verified.json'];steps.append(s)
    results,s=run_python('test_fail_closed14.py',['fail_closed14_verified.json']);negative=results['fail_closed14_verified.json'];steps.append(s)
    # Explicit cross-stage interfaces: root, local radius, global threshold,
    # full candidate graph, and the exact partition of all enumerated cases.
    rh=initial['core/root14.json'];assert root['root_sha256']==sym['root_sha256']==local['root_sha256']==forest['root_sha256']==rh
    assert forest['local_certificate_sha256']==initial['core/anchor_isolation_certificate.json']
    assert root['dimension']==116 and upper['vertices']==24 and upper['edges']==59 and upper['faces']==36
    assert upper['exact_equality_center_faces']==14 and len(upper['strict_center_faces'])==2 and len(upper['collapsed_ring_triangles'])==2
    assert F(local['anchor_Frobenius_radius'])==F(1,6) and F(local['isolation_margin'])>0
    assert partition['total']==1313024 and partition['residual']==7232 and forest['noncandidate']['cases']==7231
    excluded=sum(v['direct']+v['farkas'] for v in partition['families'].values());assert excluded==1305792
    assert excluded+forest['noncandidate']['cases']+1==partition['total']
    assert forest['candidate']['counts']=={'SPLIT':456,'DUAL':419,'LOCAL':38}
    assert forest['noncandidate']['counts']=={'DUAL':11653,'SPLIT':4422}
    assert F(forest['candidate']['minimum_force_margin'])>0 and F(forest['candidate']['minimum_local_margin'])>0 and F(forest['noncandidate']['minimum_force_margin'])>0
    final={n:sha(ROOT/n) for n in INPUTS};assert final==initial,'proof input changed during replay'
    out={'verified':True,'theorem':'r_14 = sqrt(t_star), R_14 = 1/sqrt(t_star)','global_optimality_certificate_chain_closed':True,
         'decimal_intervals':root['decimal_intervals'],'kkt_dimension':116,'root_sha256':rh,'covering_upper_verified':True,
         'enumeration_total':1313024,'metric_excluded':excluded,'noncandidate_cases':7231,'candidate_cases':1,
         'candidate_tree':forest['candidate'],'noncandidate_forest':forest['noncandidate'],'unresolved_leaves':0,
         'anchor_isolation_radius':'1/6','independent_labelled_counts':partition['labelled_totals_from_recurrence'],
         'negative_tests':negative,'input_sha256':initial,'steps':steps,'seconds':time.time()-st,
         'trust_boundary':'The written geometric/topological reductions in PROOF_zh.md, standard-library integer/Fraction arithmetic, C++17 compiler and Boost integers. No numerical optimizer is trusted. This is not a Lean/Coq formalization or a claim of external peer review.'}
    temp=MASTER.with_suffix('.tmp');temp.write_text(json.dumps(out,indent=2));os.replace(temp,MASTER)
    print('R14 COMPLETE CERTIFICATE CHAIN VERIFIED',flush=True)
    print(json.dumps({k:out[k] for k in ['verified','theorem','decimal_intervals','enumeration_total','metric_excluded','noncandidate_cases','unresolved_leaves','seconds']},indent=2),flush=True)
    return out

if __name__=='__main__':
    try:main()
    except BaseException:
        MASTER.unlink(missing_ok=True)
        raise
