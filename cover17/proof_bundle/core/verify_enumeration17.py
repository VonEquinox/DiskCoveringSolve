"""Fresh exhaustive peeling and structurally independent replay of every family.
The two implementations use different roots, region orders, fresh labels,
canonical traversals, and difference-bound closure algorithms.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import json, hashlib, tempfile, subprocess, shutil, time
import exact_force17 as G
if not __debug__: raise RuntimeError('Run without -O/PYTHONOPTIMIZE')
R=Path(__file__).resolve().parent; E=R.parent/'enumeration'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(cmd):
    p=subprocess.run(list(map(str,cmd)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    if p.returncode: raise RuntimeError('command failed: '+str(cmd)+'\n'+p.stdout)
    return p.stdout

def compile_tools(tmp):
    compiler=shutil.which('g++') or shutil.which('clang++')
    if not compiler:raise RuntimeError('A C++17 compiler is required')
    primary=tmp/'peel_enum';audit=tmp/'peel_audit'
    run([compiler,'-std=c++17','-O2','-Wall','-Wextra',E/'peel_enum.cpp','-o',primary])
    run([compiler,'-std=c++17','-O2','-Wall','-Wextra',E/'peel_audit.cpp','-o',audit])
    return primary,audit

@lru_cache(None)
def labelled_count(b,i):
    assert b>=2 and i>=0
    if b==2:return int(i==0)
    v=sum(comb(i,j)*labelled_count(k,j)*labelled_count(b-k+1,i-j)
          for k in range(2,b) for j in range(i+1))
    if i:
        v+=i*(labelled_count(b+1,i-1)-sum(comb(i-1,j)*labelled_count(3,j)*labelled_count(b,i-1-j) for j in range(i)))
    assert v>=0
    return v

def rows_from_manifest(b,path):
    rows=[]
    for idx,line in enumerate(path.read_text().splitlines()):
        ints=list(map(int,line.split()));assert len(ints)==1+32-b
        stab,*masks=ints;assert stab>0 and masks==sorted(set(masks))
        assert all(0<m<1<<17 and m.bit_count()==3 for m in masks)
        faces=sorted([[j for j in range(17) if m>>j&1] for m in masks])
        sig=G.signature(b,faces);assert sig is not None
        G.bounds(b,faces)
        rows.append({'B':b,'I':17-b,'idx':idx,'faces':faces,'signature':list(sig),'stabilizer':stab})
    return rows

def verify():
    t=time.time();results=[];rows=[];regression=[]
    with tempfile.TemporaryDirectory(prefix='cover17_enum_') as td:
        tmp=Path(td);primary,audit=compile_tools(tmp)
        for b in range(11,17):
            prefix=tmp/f'B{b}';archive=E/f'B{b}_I{17-b}.txt';report=tmp/f'audit_{b}.json'
            log=run([primary,17,b,G.A,G.C,G.C6,prefix,1])
            data=json.loads(prefix.with_suffix('.json').read_text())
            assert (data['N'],data['B'],data['I'],data['A'],data['C4'],data['C6'],data['denominator'],data['pruning'],data['complete'])==(17,b,17-b,G.A,G.C,G.C6,G.D,True,True)
            assert prefix.with_suffix('.txt').read_bytes()==archive.read_bytes(),('regenerated manifest differs',b)
            log+=run([audit,17,b,G.A,G.C,G.C6,archive,report,1])
            other=json.loads(report.read_text());assert other['verified'] and (other['N'],other['B'],other['A'],other['C4'],other['C6'])==(17,b,G.A,G.C,G.C6)
            br=rows_from_manifest(b,archive)
            assert len(br)==data['orbits']==other['surviving_orbits']
            assert data['rooted_survivors']==other['rooted_survivors']==sum(2*b//r['stabilizer'] for r in br)
            assert all(2*b%r['stabilizer']==0 for r in br)
            rows.extend(br)
            results.append({'B':b,'I':17-b,'primary':data,'independent':other,'manifest_sha256':sha(archive),'log_sha256':hashlib.sha256(log.encode()).hexdigest()})
            print('FRESH COMPLETE ENUMERATION + INDEPENDENT AUDIT',b,17-b,len(br),flush=True)
        # With all pruning disabled, compare full small labelled counts from
        # root-edge deletion, not from either implementation's traversal.
        for n in range(3,10):
            for b in range(3,n+1):
                i=n-b
                if i>3:continue
                prefix=tmp/f'small_{n}_{b}';rep=tmp/f'small_audit_{n}_{b}.json'
                run([primary,n,b,G.A,G.C,G.C6,prefix,0]);a=json.loads(prefix.with_suffix('.json').read_text())
                run([audit,n,b,G.A,G.C,G.C6,prefix.with_suffix('.txt'),rep,0]);c=json.loads(rep.read_text())
                lab=labelled_count(b,i)
                assert a['rooted_survivors']*factorial(i)==lab
                assert c['rooted_survivors']==a['rooted_survivors'] and c['surviving_orbits']==a['orbits']
                regression.append({'N':n,'B':b,'I':i,'labelled_count':lab,'rooted_count':a['rooted_survivors'],'orbits':a['orbits']})
    archived=json.loads((E/'metric17_residuals.json').read_text())
    assert rows==archived,'full residual rows do not equal regenerated surviving records'
    out={'verified':True,'N':17,'method':'exhaustive rooted peeling with safe necessary-condition pruning; independent queue/root/label/BFS/Floyd replay',
         'families':results,'surviving_orbits':len(rows),'primary_states':sum(x['primary']['calls'] for x in results),
         'independent_states':sum(x['independent']['visits'] for x in results),'unpruned_count_regressions':regression,
         'residual_sha256':sha(E/'metric17_residuals.json'),'seconds':time.time()-t}
    (R/'enumeration17_verified.json').write_text(json.dumps(out,indent=2))
    print('COMPLETE ENUMERATION VERIFIED',len(rows),'residual orbits;',len(regression),'unpruned count tests; seconds',time.time()-t,flush=True)
    return out
if __name__=='__main__':verify()
