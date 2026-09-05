"""Exact rational isolation for the 42-rod graph, in rescaled Cartesian coordinates.
No numerical root or approximate coefficient is used. A rational congruence
is only a proposal: all of its acceptance checks are exact integer operations.
"""
from fractions import Fraction as F
from pathlib import Path
from math import lcm
import json,time,hashlib
import geometry19 as S
if not __debug__:raise RuntimeError('Run without -O')
R=S.R; GAMMA=400; LAMBDA=F(1,650); NU=F(1,156); RADIUS=F(1,5);QR=10**10

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
 S.check_geometry();A,M,P=S.matrices();g=S.G
 di=[w*w for w in S.WEIGHTS]+[NU*NU]*(S.B-1)
 assert all(-NU==m for m in S.MUS)
 H=[row.copy() for row in M]
 for k,row in enumerate(A):
  nz=[(i,v) for i,v in enumerate(row) if v]
  for i,v in nz:
   for j,w in nz:H[i][j]+=GAMMA*di[k]*v*w
 for i in range(g):H[i][i]-=LAMBDA*P[i]
 den=lcm(*(v.denominator for row in H for v in row))
 ints=[[(v*den).numerator for v in row] for row in H]
 assert LAMBDA-2*GAMMA*NU*NU*RADIUS*RADIUS>0
 return ints,den

def verify():
 st=time.time();H,den=prepare();g=S.G
 path=R/'anchor_isolation_certificate.json';cert=json.load(open(path))
 assert cert['QR']==QR and cert['geometry_sha256']==sha(R/'geometry19.py')
 T=cert['Rnum'];assert len(T)==g and all(len(row)==g for row in T)
 assert all(type(v)is int for row in T for v in row)
 assert all(T[i][j]==0 for i in range(g) for j in range(i)) and all(T[i][i]!=0 for i in range(g))
 HT=[[sum(H[i][k]*T[k][j] for k in range(j+1)) for j in range(g)] for i in range(g)]
 C=[[sum(T[k][i]*HT[k][j] for k in range(i+1)) for j in range(g)] for i in range(g)]
 assert all(C[i][j]==C[j][i] for i in range(g) for j in range(g))
 margins=[C[i][i]-sum(abs(C[i][j]) for j in range(g) if i!=j) for i in range(g)]
 assert min(margins)>0
 out={'verified':True,'dimension':g,'rods':S.M,'gamma':GAMMA,'lambda':str(LAMBDA),'nu':str(NU),'anchor_radius':str(RADIUS),'isolation_margin':str(LAMBDA-2*GAMMA*NU**2*RADIUS**2),'congruence_margin':str(F(min(margins),den*QR**2)),'geometry_sha256':sha(R/'geometry19.py'),'certificate_sha256':sha(path),'seconds':time.time()-st,'conclusion':'All nodes of the 42-rod subgraph equal the candidate; center c12 and unused face nodes are not constrained by this theorem.'}
 (R/'anchor_isolation_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return out
if __name__=='__main__':verify()
