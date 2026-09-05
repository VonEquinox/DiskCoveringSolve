#!/usr/bin/env python3
"""Full replay. No numerical optimizer or pre-existing success report is trusted."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Proof verification requires assertions; do not use -O/-OO/PYTHONOPTIMIZE')
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent
CORE=ROOT/'core'; ENUM=ROOT/'enumeration'; WORK=ROOT/'replay_workspace'; REP=ROOT/'reports'
FAMILIES=[(b,15-b) for b in range(10,15)]
MASTER=REP/'MASTER_VERIFIED.json'

def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def input_hashes():
    paths=[ROOT/'verify_all.py',ROOT/'PROOF_zh.md',ROOT/'README.md']
    paths += [p for folder in (CORE,ENUM) for p in folder.iterdir()
              if p.is_file() and not p.name.endswith('_verified.json') and
              (p.suffix in ('.py','.cpp','.json','.gz'))]
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--regenerate',action='store_true',help='Regenerate all representative files before the independent census audit')
    ap.add_argument('--jobs',type=int,default=min(4,os.cpu_count() or 1),help='Maximum concurrent C++ generation/audit processes')
    args=ap.parse_args()
    if not 1<=args.jobs<=16:ap.error('--jobs must lie between 1 and 16')
    WORK.mkdir(exist_ok=True);REP.mkdir(exist_ok=True);MASTER.unlink(missing_ok=True)
    started=time.time();hash_before=input_hashes()
    for name in ['root15_verified.json','upper15_verified.json','anchor_isolation_verified.json','partition15_verified.json','forest15_verified.json','rejection15_verified.json']:
        (CORE/name).unlink(missing_ok=True)
    for B,_ in FAMILIES:(REP/f'audit{B}.json').unlink(missing_ok=True)
    full=(REP/'FULL_REPLAY.log').open('w',encoding='utf8')
    def note(text):
        print(text,flush=True);print(text,file=full,flush=True)
    def run(label,cmd,cwd=ROOT):
        note('RUN '+label)
        lp=REP/(label+'.log')
        with lp.open('w',encoding='utf8') as log:
            rr=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
        if rr.returncode:
            note(lp.read_text(encoding='utf8',errors='replace')[-10000:])
            raise RuntimeError(f'{label}: exit {rr.returncode}; see {lp}')
        txt=lp.read_text(encoding='utf8',errors='replace')
        print(txt,file=full,end='' if txt.endswith('\n') else '\n',flush=True)
        note('PASS '+label)
    try:
        compiler=os.environ.get('CXX') or shutil.which('g++') or shutil.which('clang++')
        if not compiler:raise RuntimeError('A C++17 compiler (g++ or clang++) is required')
        for name in ['enumerate15','audit15','metric_audit15']:
            run('compile_'+name,[compiler,'-O3','-std=c++17',ENUM/(name+'.cpp'),'-o',WORK/name])
        # Generation is untrusted. Every cached or newly generated file is audited below.
        def gen(fam):
            B,I=fam;prefix=WORK/f'B{B}_I{I}'
            if args.regenerate or not prefix.with_suffix('.txt').exists():
                run(f'generate_B{B}',[WORK/'enumerate15',B,I,prefix])
            return B
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for fut in as_completed([pool.submit(gen,f) for f in FAMILIES]):fut.result()
        def audit(fam):
            B,I=fam
            run(f'audit_B{B}',[WORK/'audit15',WORK/f'B{B}_I{I}.txt',REP/f'audit{B}.json'])
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for fut in as_completed([pool.submit(audit,f) for f in FAMILIES]):fut.result()
        for name in ['verify_root15.py','verify_upper15.py','anchor_isolation15.py','verify_partition15.py','verify_forest15.py','test_fail_closed15.py']:
            run(name[:-3],[sys.executable,'-S','-B',CORE/name])
        def load(name):return json.loads((CORE/name).read_text())
        root=load('root15_verified.json');upper=load('upper15_verified.json');local=load('anchor_isolation_verified.json');part=load('partition15_verified.json');forest=load('forest15_verified.json');reject=load('rejection15_verified.json')
        assert all(r['verified'] is True for r in [root,upper,local,part,forest,reject])
        assert root['dimension']==148 and root['root_sha256']==sha(CORE/'root15.json')
        assert local['root_sha256']==root['root_sha256']==forest['root_sha256']
        assert local['certificate_sha256']==forest['local_certificate_sha256']==sha(CORE/'anchor_isolation_certificate.json')
        assert (upper['vertices'],upper['edges'],upper['faces'],upper['boundary_caps'],upper['exact_equality_center_faces'])==(26,64,39,11,13)
        assert len(upper['strict_center_faces'])==4
        assert part['total']==11950884 and part['excluded']==11907870 and part['residual']==43014
        nc=forest['noncandidate'];cand=forest['candidate'];cc=cand['counts'];nn=nc['counts']
        assert nc['cases']==43013 and part['excluded']+nc['cases']+1==part['total']
        assert cc=={'SPLIT':491,'DUAL':459,'LOCAL':33} and cand['nodes']==983
        assert set(nn)<= {'SPLIT','DUAL','EMPTY'} and sum(nn.values())==nc['nodes']
        assert nn.get('DUAL',0)+nn.get('EMPTY',0)==nn.get('SPLIT',0)+nc['cases']
        assert all(nc['families'][str(B)]['cases']==part['families'][f'B{B}_I{I}']['residual']-(B==11) for B,I in FAMILIES)
        assert input_hashes()==hash_before,'Proof input changed during replay'
        for name in ['root15_verified.json','upper15_verified.json','anchor_isolation_verified.json','partition15_verified.json','forest15_verified.json','rejection15_verified.json']:
            shutil.copy2(CORE/name,REP/name)
        out={'version':1,'verified':True,'global_optimality_certificate_chain_closed':True,
             'statement':'r_15 = sqrt(t_star), where t_star is the t-coordinate of the unique root of core/system15.py in core/root15.json',
             'root_dimension':148,'decimal_intervals':root['decimal_intervals'],
             'enumeration_total':part['total'],'metric_excluded':part['excluded'],
             'noncandidate_cases':nc['cases'],'candidate_cases':1,'unresolved_leaves':0,
             'candidate_tree':cand,'noncandidate_forest':nc,
             'upper_cover':{k:upper[k] for k in ['vertices','edges','faces','boundary_caps','exact_equality_center_faces']},
             'anchor_isolation':{k:local[k] for k in ['gamma','lambda','nu_upper','anchor_Frobenius_radius','isolation_margin']},
             'enumeration_audits':{str(B):json.loads((REP/f'audit{B}.json').read_text()) for B,_ in FAMILIES},
             'rejection_tests':reject,'proof_input_sha256':hash_before,
             'completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.time()-started,
             'external_peer_review_completed':False,'proof_assistant_formalization':False}
        MASTER.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
        (ROOT/'SHA256SUMS_PROOF_INPUTS.txt').write_text(''.join(h+'  '+p+'\n' for p,h in hash_before.items()),encoding='utf8')
        note('R15 COMPLETE CERTIFICATE CHAIN VERIFIED')
        note(f'11907870 + 43013 + 1 = 11950884; unresolved leaves = 0')
        return out
    except BaseException:
        MASTER.unlink(missing_ok=True)
        note('R15 VERIFICATION FAILED; no master success report retained')
        raise
    finally:full.close()

if __name__=='__main__':main()
