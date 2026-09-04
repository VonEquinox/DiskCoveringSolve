from collections import defaultdict,deque
import itertools,time,pickle,sys

def canon_face(f):return tuple(sorted(f))
def canon_state(fs):return tuple(sorted(canon_face(f) for f in fs))

def polygon_fan(B):return canon_state((0,i,i+1) for i in range(1,B-1))
def insert_vertex(state,v,face_index=0):
    fs=list(state);a,b,c=fs.pop(face_index)
    fs += [(v,a,b),(v,b,c),(v,c,a)]
    return canon_state(fs)

def flips(state,B):
    em=defaultdict(list)
    edges=set()
    for fi,f in enumerate(state):
        a,b,c=f
        for x,y in ((a,b),(a,c),(b,c)):
            e=(x,y) if x<y else (y,x);em[e].append(fi);edges.add(e)
    for (a,b),idxs in em.items():
        if len(idxs)!=2:continue
        f1=state[idxs[0]];f2=state[idxs[1]]
        c=next(x for x in f1 if x not in (a,b));d=next(x for x in f2 if x not in (a,b))
        if c==d:continue
        cd=(c,d) if c<d else (d,c)
        if cd in edges:continue
        fs=[f for k,f in enumerate(state) if k not in idxs]
        fs += [(c,d,a),(c,d,b)]
        ns=canon_state(fs)
        # basic manifold checks
        yield ns

def enumerate_all(B,I,maxstates=2000000):
    st=polygon_fan(B)
    for j in range(I):
        st=insert_vertex(st,B+j,0)
    seen={st};q=deque([st]);t=time.time();nflip=0
    while q:
        s=q.popleft()
        for z in flips(s,B):
            nflip+=1
            if z not in seen:
                seen.add(z);q.append(z)
                if len(seen)%10000==0:print('states',len(seen),'queue',len(q),'sec',time.time()-t,flush=True)
                if len(seen)>maxstates:raise RuntimeError('too many')
    print('DONE B I',B,I,'states',len(seen),'attempts',nflip,'sec',time.time()-t)
    return seen

def transform_state(s,B,I,shift=0,reflect=False,swap=False):
    def mp(x):
        if x<B:
            y=(-x if reflect else x)+shift
            return y%B
        if I==2 and swap:return B+(1-(x-B))
        return x
    return canon_state(tuple(mp(x) for x in f) for f in s)

def canonical_orbit(s,B,I):
    arr=[]
    for refl in (False,True):
      for sh in range(B):
       arr.append(transform_state(s,B,I,sh,refl,False))
       if I==2:arr.append(transform_state(s,B,I,sh,refl,True))
    return min(arr)

def main(B,I):
    states=enumerate_all(B,I)
    unseen=set(states);reps=[]
    while unseen:
        s=next(iter(unseen));orb=set()
        for refl in (False,True):
          for sh in range(B):
            orb.add(transform_state(s,B,I,sh,refl,False))
            if I==2:orb.add(transform_state(s,B,I,sh,refl,True))
        unseen.difference_update(orb);reps.append(min(orb))
    reps=set(reps)
    print('orbits',len(reps))
    out=f'/mnt/data/triangulations_B{B}_I{I}.pkl'
    with open(out,'wb') as f:pickle.dump(sorted(reps),f)
    print(out)
if __name__=='__main__':main(int(sys.argv[1]),int(sys.argv[2]))
