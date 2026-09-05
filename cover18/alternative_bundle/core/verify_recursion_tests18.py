"""Unpruned small-case crosscheck against an independently derived recurrence.

These tests are additional checks, not the completeness proof of the n=18 run.
"""
from pathlib import Path
from math import factorial
import json, subprocess, tempfile, time, hashlib
from verify_peeling18 import labelled
if not __debug__:
    raise RuntimeError('Assertions must be enabled')
ROOT=Path(__file__).resolve().parents[1]

def verify():
    start=time.time(); tests=[]
    cases=[(n,b) for n in range(5,11) for b in range(max(3,n-4),n+1)]
    cases += [(11,8),(12,9)]
    with tempfile.TemporaryDirectory(prefix='recursion18_',dir=ROOT/'_runtime') as tmp:
        for n,b in cases:
            outputs=[]
            for method in ('dfs','bfs'):
                path=Path(tmp)/f'{method}_{n}_{b}.json'
                args=[str(ROOT/'_runtime'/f'peel_{method}18'),str(n),str(b),str(path)]
                args+=['0','noscreen'] if method=='dfs' else ['no-screen']
                p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
                assert p.returncode==0,(args,p.stdout,p.stderr)
                d=json.loads(path.read_text());assert d['screen'] is False
                assert d['N']==n and d['B']==b and d['I']==n-b
                outputs.append(d)
            a,c=outputs; expected=labelled(b,n-b)//factorial(n-b)
            assert a['complete_rooted']==c['rooted_survivors']==expected
            assert a['distinct_survivors']==c['surviving_orbits']
            tests.append({'N':n,'B':b,'I':n-b,'labelled':labelled(b,n-b),
                          'rooted':expected,'orbits':c['surviving_orbits']})
    out={'verified':True,'tests':tests,'cases':len(tests),'pruning_disabled':True,
         'role':'Additional finite checks; n=18 completeness is the root-face induction plus freshly executed exhaustive searches.',
         'seconds':time.time()-start}
    (ROOT/'reports'/'recursion_tests18_verified.json').write_text(json.dumps(out,indent=2))
    print('UNPRUNED INDEPENDENT RECURRENCE CHECKS VERIFIED',len(tests),'families',flush=True)
    return out
if __name__=='__main__':verify()
