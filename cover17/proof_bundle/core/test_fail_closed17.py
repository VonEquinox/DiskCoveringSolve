"""Adversarial regression tests. These supplement, not replace, the proof."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from contextlib import contextmanager, redirect_stdout
import copy, gzip, json, io, subprocess, tempfile, shutil, sys, time
import system17 as S
import exact_force17 as G
import verify_forest17 as V
import verify_root17 as RV
import anchor_isolation17 as L
import verify_interfaces17 as IF
from verify_enumeration17 import compile_tools, run
if not __debug__:raise RuntimeError('Run without -O/PYTHONOPTIMIZE')
R=Path(__file__).resolve().parent;E=R.parent/'enumeration'

@contextmanager
def patch(obj, name, value):
    old=getattr(obj,name);setattr(obj,name,value)
    try:yield
    finally:setattr(obj,name,old)

def main():
    t=time.time();accepted=[]
    def reject(name,fn,allow_nonpositive=False):
        try:
            with redirect_stdout(io.StringIO()):result=fn()
        except (AssertionError,RuntimeError,ValueError,KeyError,IndexError,TypeError,FileNotFoundError,OverflowError):
            accepted.append(name);return
        if allow_nonpositive and (result is None or (isinstance(result,F) and result<=0)):
            accepted.append(name);return
        raise AssertionError('corruption was accepted: '+name)
    def must_fail_cmd(name,args):
        p=subprocess.run(list(map(str,args)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        assert p.returncode!=0,(name,p.stdout)
        accepted.append(name)
    with gzip.open(E/'noncandidate17_trees.jsonl.gz','rt') as f:
        head=json.loads(next(f));first=json.loads(next(f))
    assert len(first['nodes'])==1 and first['nodes'][0]['kind']=='DUAL'
    b=first['B'];fa=first['faces'];lo,hi=G.rootbox(G.bounds(b,fa));lo,hi=G.tighten(G.bounds(b,fa),lo,hi);ed=G.edges(b,fa)
    leaf=first['nodes'][0]
    assert G.check(b,ed,lo,hi,leaf)>0
    bad=copy.deepcopy(leaf);bad['force_num'][0][0]+=1
    reject('one integer rod-force unit changed',lambda:G.check(b,ed,lo,hi,bad))
    bad=copy.deepcopy(leaf);bad['lambda_num'][0]=-1
    reject('negative arc-halfplane multiplier',lambda:G.check(b,ed,lo,hi,bad))
    bad=copy.deepcopy(leaf);bad['force_num'].pop()
    reject('missing rod force',lambda:G.check(b,ed,lo,hi,bad))
    bad=copy.deepcopy(leaf);bad['lambda_num'][0]=float(bad['lambda_num'][0])
    reject('floating-point integer certificate field',lambda:G.check(b,ed,lo,hi,bad))
    bad=copy.deepcopy(leaf);bad['force_num']=[[0,0] for _ in ed];bad['support_num']=[[0,0] for _ in range(b-1)];bad['lambda_num']=[0]*(b-1)
    reject('zero total force gives no strict separation',lambda:G.check(b,ed,lo,hi,bad),True)
    reject('unresolved leaf',lambda:V.traverse(b,fa,[{'kind':'UNRESOLVED'}],False))
    reject('nonempty box falsely marked EMPTY',lambda:V.traverse(b,fa,[{'kind':'EMPTY'}],False))
    reject('LOCAL leaf in a noncandidate tree',lambda:V.traverse(b,fa,[{'kind':'LOCAL'}],False))
    reject('whole candidate root box falsely marked LOCAL',lambda:V.traverse(S.B,S.FACES,[{'kind':'LOCAL'}],True))
    # A one-node valid tree with appended garbage must fail reachability.
    reject('unreachable appended node',lambda:V.traverse(b,fa,[copy.deepcopy(leaf),{'kind':'EMPTY'}],False))
    boxlo,boxhi=G.rootbox(G.bounds(b,fa));boxlo,boxhi=G.tighten(G.bounds(b,fa),boxlo,boxhi)
    j=next(i for i,(a,c) in enumerate(zip(boxlo,boxhi)) if a<c);mid=str((boxlo[j]+boxhi[j])/2)
    reject('duplicate split children',lambda:V.traverse(b,fa,[{'kind':'SPLIT','axis':j,'mid':mid,'children':[1,1]},leaf],False))
    reject('cyclic split edge',lambda:V.traverse(b,fa,[{'kind':'SPLIT','axis':j,'mid':mid,'children':[0,1]},leaf],False))
    reject('split at the boundary instead of the interior',lambda:V.traverse(b,fa,[{'kind':'SPLIT','axis':j,'mid':str(boxlo[j]),'children':[1,2]},leaf,leaf],False))
    rows=json.loads((E/'metric17_residuals.json').read_text())
    with tempfile.TemporaryDirectory(prefix='cover17_negative_') as td:
        tmp=Path(td);enum=tmp/'enumeration';core=tmp/'core';enum.mkdir();core.mkdir()
        # The map verifier must not accept an arbitrary permutation or record number.
        mp=json.loads((E/'candidate17_map.json').read_text());p=copy.deepcopy(mp);p['map'][0],p['map'][1]=p['map'][1],p['map'][0]
        (enum/'candidate17_map.json').write_text(json.dumps(p))
        with patch(V,'ENUM',enum):reject('incorrect candidate graph isomorphism',lambda:V.check_mapping(rows))
        (enum/'candidate17_map.json').write_text(json.dumps(mp))
        with patch(V,'ENUM',enum):reject('missing candidate residual record',lambda:V.check_mapping(rows[:mp['residual_index']]))
        # Corrupt the root in a private temporary directory, never the proof inputs.
        root=json.loads((R/'root17.json').read_text());r=copy.deepcopy(root);r['xnum'][0]=str(int(r['xnum'][0])+10**85)
        (core/'root17.json').write_text(json.dumps(r))
        with patch(S,'R',core),patch(RV,'R',core):reject('algebraic-root midpoint displaced',RV.verify)
        (core/'root17.json').write_text(json.dumps(root));shutil.copy2(R/'system17.py',core/'system17.py')
        # Preserve original bytes so the authentic root hash remains valid.
        shutil.copy2(R/'root17.json',core/'root17.json')
        cert=json.loads((R/'anchor_isolation_certificate.json').read_text());cert['Rnum'][0][0]=0
        (core/'anchor_isolation_certificate.json').write_text(json.dumps(cert))
        with patch(L,'ROOT',core):reject('singular alleged positive-definite congruence',L.main)
        with patch(G,'C6',35300):reject('unverified six-rod angle cap',IF.verify)
        with patch(L,'RADIUS',F(1,2)):reject('unlinked inflated local radius',IF.verify)
        # Missing and duplicated noncandidate records use the same production reader.
        for f in ('metric17_residuals.json','candidate17_map.json'):shutil.copy2(E/f,enum/f)
        shutil.copy2(R/'candidate17_tree.json.gz',core/'candidate17_tree.json.gz')
        shutil.copy2(R/'anchor_isolation_certificate.json',core/'anchor_isolation_certificate.json')
        with gzip.open(enum/'noncandidate17_trees.jsonl.gz','wt') as f:f.write(json.dumps(head)+'\n')
        with patch(V,'R',core),patch(V,'ENUM',enum):reject('missing all noncandidate records despite correct header',V.verify)
        with gzip.open(enum/'noncandidate17_trees.jsonl.gz','wt') as f:
            f.write(json.dumps(head)+'\n'+json.dumps(first)+'\n'+json.dumps(first)+'\n')
        with patch(V,'R',core),patch(V,'ENUM',enum):reject('duplicate noncandidate record',V.verify)
        # Independently rebuilt C++ auditor: invalid/missing/duplicate manifests.
        primary,audit=compile_tools(tmp);manifest=(E/'B16_I1.txt').read_text();args=[audit,17,16,G.A,G.C,G.C6,tmp/'bad.txt',tmp/'bad.json',1]
        def cppbad(name,text):
            (tmp/'bad.txt').write_text(text);(tmp/'bad.json').unlink(missing_ok=True)
            must_fail_cmd(name,args);assert not (tmp/'bad.json').exists()
        lines=manifest.splitlines()
        cppbad('duplicate valid surviving orbit',manifest+lines[0]+'\n')
        cppbad('omitted surviving orbit', '\n'.join(lines[1:])+'\n')
        cppbad('invalid triangle mask',lines[0]+' 1\n')
        cppbad('trailing nonnumeric manifest bytes',manifest+'garbage\n')
        cppbad('overflowed integer manifest token',manifest+lines[0]+' 18446744073709551616\n')
        (tmp/'bad.txt').unlink();must_fail_cmd('missing enumeration manifest',args)
        prefix=tmp/'limited';prefix.with_suffix('.json').write_text('{"complete":true}');prefix.with_suffix('.txt').write_text(manifest)
        must_fail_cmd('interrupted peeling cannot leave stale success',[primary,17,11,G.A,G.C,G.C6,prefix,1,1])
        assert not prefix.with_suffix('.json').exists() and not prefix.with_suffix('.txt').exists()
        # Master guard is tested in isolation; it cannot erase the actual replay output.
        shutil.copy2(R.parent/'verify_all.py',tmp/'verify_all.py');(tmp/'reports').mkdir();(tmp/'reports/MASTER_VERIFIED.json').write_text('{"verified":true}')
        must_fail_cmd('disabled assertions in the master',[sys.executable,'-O','-S','-B',tmp/'verify_all.py'])
        assert not (tmp/'reports/MASTER_VERIFIED.json').exists()
    must_fail_cmd('disabled assertions in the force checker',[sys.executable,'-O','-S','-B',R/'exact_force17.py'])
    out={'verified':True,'tests_passed':len(accepted),'tests':accepted,'seconds':time.time()-t}
    (R/'negative17_verified.json').write_text(json.dumps(out,indent=2))
    print('FAIL-CLOSED TESTS VERIFIED',len(accepted),'seconds',time.time()-t,flush=True)
    return out
if __name__=='__main__':main()
