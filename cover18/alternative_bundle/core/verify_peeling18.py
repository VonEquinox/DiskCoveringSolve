"""Link exhaustive partial-topology searches, exact angular bounds, and survivors.

The two C++ programs are re-executed by verify_all.py. They do not merely read
an exclusion manifest. This module checks both outputs against the very same
complete graphs used by all angle-tree certificates.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import time
import system18 as S
import exact_force18 as G
from exact_arcs import sincos_turn, SCALE
if not __debug__:
    raise RuntimeError('Exact verification requires assertions enabled')
ROOT=Path(__file__).resolve().parents[1]
ENUM=ROOT/'enumeration'
RUNTIME=ROOT/'_runtime'
REPORTS=ROOT/'reports'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

@lru_cache(None)
def labelled(B,I):
    """Independent root-edge recurrence, with labelled interior vertices."""
    assert type(B) is int and type(I) is int and B>=2 and I>=0
    if B==2:
        return int(I==0)
    ans=sum(comb(I,j)*labelled(k,j)*labelled(B-k+1,I-j)
            for k in range(2,B) for j in range(I+1))
    if I:
        ans+=I*(labelled(B+1,I-1)-sum(
            comb(I-1,j)*labelled(3,j)*labelled(B,I-1-j)
            for j in range(I)))
    return ans

def check_constants():
    assert S.N==G.N==18 and S.B==12 and G.V==3*S.N-2
    assert (G.A,G.C,G.C6,G.D)==(937111,1970786,3362065,10000000)
    d=json.loads((S.R/'root18.json').read_text())
    t=F(int(d['xnum'][S.TID]),int(d['Qx']));rho=F(1,10**80)
    assert 0<t-rho<t+rho<G.RU**2<F(1,4)
    gaps=[F(G.A,G.D),F(G.C,G.D),F(G.C6,G.D)]
    assert 0<gaps[0]<gaps[1]<gaps[2]<F(1,2)
    assert 10*gaps[0]<1
    margins=[]
    for k,cap in enumerate(gaps,1):
        _,sin_iv=sincos_turn(cap/2)
        margin=F(sin_iv[0],SCALE)-k*G.RU
        assert margin>0
        margins.append(str(margin))
    return {'radius_upper':str(G.RU),'turn_caps':list(map(str,gaps)),
            'sine_margins':margins,'minimum_boundary_cells':11,
            'root_sha256':sha(S.R/'root18.json')}

def bfs_code(B,faces):
    """Independent Python code: inspect all oriented roots, not a hash guess."""
    N=18
    assert 11<=B<N and len(faces)==2*N-B-2
    faces=[tuple(f) for f in faces]
    assert len(set(faces))==len(faces)
    assert all(len(f)==3 and list(f)==sorted(set(f)) and 0<=f[0]<f[-1]<N for f in faces)
    inc={}
    for k,f in enumerate(faces):
        for i in range(3):
            e=tuple(sorted((f[i],f[(i+1)%3])))
            inc.setdefault(e,[]).append(k)
    bd={tuple(sorted((i,(i+1)%B))) for i in range(B)}
    assert all(len(fs)==(1 if e in bd else 2) for e,fs in inc.items())
    assert bd<=inc.keys()
    best=None;stab=0
    for sign in (1,-1):
        for start in range(B):
            image={(start+sign*i)%B:i for i in range(B)}
            u,v=start,(start+sign)%B
            root=inc[tuple(sorted((u,v)))];assert len(root)==1
            queue=deque([(root[0],u,v)]);scheduled={root[0]};code=[]
            while queue:
                f,a,b=queue.popleft()
                w=next(x for x in faces[f] if x!=a and x!=b)
                if w not in image:
                    image[w]=len(image)
                assert a in image and b in image
                code.append((1<<image[a])|(1<<image[b])|(1<<image[w]))
                for edge in ((b,w),(w,a)):
                    for j in inc[tuple(sorted(edge))]:
                        if j not in scheduled:
                            scheduled.add(j);queue.append((j,*edge))
            assert len(image)==N and len(code)==len(faces)
            code=tuple(sorted(code))
            if best is None or code<best:
                best=code;stab=1
            elif code==best:
                stab+=1
    return best,stab

def verify():
    started=time.time();constants=check_constants()
    rows=json.loads((ENUM/'survivors18.json').read_text())
    assert isinstance(rows,list) and rows
    cursor=0;families={};total=0;rooted=0
    for B in range(11,18):
        primary=json.loads((RUNTIME/f'dfs_B{B}.json').read_text())
        audit=json.loads((RUNTIME/f'bfs_B{B}.json').read_text())
        for out in (primary,audit):
            assert out['N']==18 and out['B']==B and out['I']==18-B
            assert out['screen'] is True
            assert out['angle_scale']==G.D and out['angle_caps']==[G.A,G.C,G.C6]
        assert audit['verified'] is True
        representatives=primary['survivors']
        assert len(representatives)==primary['distinct_survivors']==audit['surviving_orbits']
        assert primary['complete_rooted']==audit['rooted_survivors']
        expected={}
        for index,masks in enumerate(representatives):
            assert all(type(m) is int and 0<m<(1<<18) and m.bit_count()==3 for m in masks)
            assert masks==sorted(set(masks))
            fs=[[i for i in range(18) if m>>i&1] for m in masks]
            row={'B':B,'I':18-B,'idx':index,'faces':fs}
            assert cursor<len(rows) and row==rows[cursor]
            # Rebuild the necessary angular domain from the final face set.
            mat=G.bounds(B,fs)
            assert len(mat)==B and all(mat[i][i]==0 for i in range(B))
            code,stabilizer=bfs_code(B,fs)
            assert code not in expected
            expected[code]=stabilizer
            cursor+=1
        observed={}
        for rec in audit['orbits']:
            key=tuple(rec['masks']);stab=rec['stabilizer']
            assert key not in observed and type(stab) is int and stab>0 and (2*B)%stab==0
            assert rec['rooted_multiplicity']==2*B//stab
            observed[key]=stab
        assert observed==expected
        assert sum(2*B//s for s in expected.values())==audit['rooted_survivors']
        count=labelled(B,18-B)
        assert count%factorial(18-B)==0
        closed=2*factorial(2*B-3)*factorial(4*(18-B)+2*B-5)//(
            factorial(B-1)*factorial(B-3)*factorial(3*(18-B)+2*B-3))
        assert count==closed
        families[f'B{B}_I{18-B}']={
            'dfs_calls':primary['calls'],'dfs_metric_cuts':primary['prunes'],
            'bfs_calls':audit['calls'],'bfs_metric_cuts':audit['metric_rejections'],
            'rooted_survivors':audit['rooted_survivors'],'surviving_orbits':len(expected),
            'labelled_unpruned_count_from_recurrence':count,
            'unlabelled_interior_rooted_count':count//factorial(18-B),
            'dfs_output_sha256':sha(RUNTIME/f'dfs_B{B}.json'),
            'bfs_output_sha256':sha(RUNTIME/f'bfs_B{B}.json')}
        total+=len(expected);rooted+=audit['rooted_survivors']
        print('PARTIAL TOPOLOGY FAMILY VERIFIED',B,len(expected),audit['rooted_survivors'],flush=True)
    assert cursor==len(rows)==total
    out={'verified':True,'method':'Exhaustive root-face recursion with exact cuts BEFORE completion; independent DFS/Bellman and BFS/Floyd implementations.',
         'constants':constants,'families':families,'surviving_orbits':total,'rooted_survivors':rooted,
         'all_families':[11,12,13,14,15,16,17],'unresolved_combinatorial_branches':0,
         'survivors_sha256':sha(ENUM/'survivors18.json'),'seconds':time.time()-started}
    (REPORTS/'peeling18_verified.json').write_text(json.dumps(out,indent=2))
    print('EXHAUSTIVE PARTIAL TOPOLOGY PARTITION VERIFIED',total,'orbits',flush=True)
    return out
if __name__=='__main__':
    verify()
