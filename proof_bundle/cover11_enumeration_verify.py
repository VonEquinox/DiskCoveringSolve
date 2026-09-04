#!/usr/bin/env python3
import sys,pickle,time,hashlib,json,math
sys.path.insert(0,'/mnt/data')
import enumerate_disk_triangulations as E

def brown_labeled_count(B,I):
 m=B-3;n=I
 rooted=(2*math.factorial(2*m+3)*math.factorial(4*n+2*m+1)//
         (math.factorial(m+2)*math.factorial(m)*math.factorial(n)*math.factorial(3*n+2*m+3)))
 return rooted*math.factorial(I)

def verify(B,I):
 t=time.time();states=E.enumerate_all(B,I)
 brown=brown_labeled_count(B,I);assert len(states)==brown,(B,I,len(states),brown)
 reps={E.canonical_orbit(s,B,I) for s in states}
 stored=pickle.load(open(f'/mnt/data/triangulations_B{B}_I{I}.pkl','rb'))
 assert sorted(reps)==stored,(B,I,len(reps),len(stored))
 raw=pickle.dumps(stored,protocol=4)
 out={'B':B,'I':I,'labeled_states':len(states),'brown_labeled_count':brown,'brown_equal':True,'orbits':len(reps),'stored_equal':True,'sha256_protocol4':hashlib.sha256(raw).hexdigest(),'sec':time.time()-t}
 print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 xs=[verify(9,1),verify(9,2),verify(10,1)]
 assert [x['orbits'] for x in xs]==[291,2548,1004]
 open('/mnt/data/cover11_enumeration_verify.json','w').write(json.dumps(xs,indent=2))
