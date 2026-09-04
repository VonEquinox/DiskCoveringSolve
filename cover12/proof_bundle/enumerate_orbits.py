#!/usr/bin/env python3
from collections import defaultdict, deque
import itertools, pickle, sys, time, math
from pathlib import Path
from fractions import Fraction

def canon_face(f): return tuple(sorted(f))
def canon_state(fs): return tuple(sorted(canon_face(f) for f in fs))
def polygon_fan(B): return canon_state((0,i,i+1) for i in range(1,B-1))
def insert_vertex(state,v,face_index=0):
    fs=list(state); a,b,c=fs.pop(face_index); fs += [(v,a,b),(v,b,c),(v,c,a)]; return canon_state(fs)
def flips(state):
    em=defaultdict(list); edges=set()
    for fi,f in enumerate(state):
        a,b,c=f
        for x,y in ((a,b),(a,c),(b,c)):
            e=(x,y) if x<y else (y,x); em[e].append(fi); edges.add(e)
    for (a,b),idxs in em.items():
        if len(idxs)!=2: continue
        f1,f2=state[idxs[0]],state[idxs[1]]
        c=next(x for x in f1 if x not in (a,b)); d=next(x for x in f2 if x not in (a,b))
        if c==d: continue
        cd=(c,d) if c<d else (d,c)
        if cd in edges: continue
        fs=[f for k,f in enumerate(state) if k not in idxs]
        fs += [(c,d,a),(c,d,b)]
        yield canon_state(fs)
def transform_state(s,B,I,shift,reflect,perm):
    def mp(x):
        if x<B: return ((-x if reflect else x)+shift)%B
        return B+perm[x-B]
    return canon_state(tuple(mp(x) for x in f) for f in s)
def canonical_orbit(s,B,I):
    best=None; perms=list(itertools.permutations(range(I))) if I else [()]
    for reflect in (False,True):
        for shift in range(B):
            for perm in perms:
                z=transform_state(s,B,I,shift,reflect,perm)
                if best is None or z<best: best=z
    return best
def enumerate_orbits(B,I,maxorbits=1000000):
    st=polygon_fan(B)
    for j in range(I): st=insert_vertex(st,B+j,0)
    st=canonical_orbit(st,B,I); seen={st}; q=deque([st]); t=time.time(); attempts=0
    while q:
        s=q.popleft()
        for z in flips(s):
            attempts += 1; z=canonical_orbit(z,B,I)
            if z not in seen:
                seen.add(z); q.append(z)
                if len(seen)%1000==0:
                    print('B,I',B,I,'orbits',len(seen),'queue',len(q),'attempts',attempts,'sec',time.time()-t,flush=True)
                if len(seen)>maxorbits: raise RuntimeError('too many orbits')
    print('DONE B,I',B,I,'orbits',len(seen),'attempts',attempts,'sec',time.time()-t,flush=True)
    return sorted(seen)
def brown_labeled(B,I):
    m=B-3; n=I
    x=Fraction(2*math.factorial(2*m+3)*math.factorial(4*n+2*m+1),math.factorial(m+2)*math.factorial(m)*math.factorial(n)*math.factorial(3*n+2*m+3))
    return int(x*math.factorial(n))
if __name__=='__main__':
    B,I=map(int,sys.argv[1:3]); out=enumerate_orbits(B,I)
    fn=str(Path(__file__).resolve().parent/f'triangulations_B{B}_I{I}_orbits.pkl'); pickle.dump(out,open(fn,'wb'),protocol=5)
    print('saved',fn,'brown_labeled',brown_labeled(B,I))
