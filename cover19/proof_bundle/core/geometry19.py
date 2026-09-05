"""Exact triangular-lattice candidate in Q(sqrt(3)); standard library only.
A point (x,y) in this file means Cartesian (x, sqrt(3)*y).
"""
from fractions import Fraction as F
from itertools import combinations
B=12; N=19; T=F(1,13)
AB=[(2,0),(1,1),(0,2),(-1,2),(-2,2),(-2,1),(-2,0),(-1,-1),(0,-2),(1,-2),(2,-2),(2,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0)]
RAW=[(F(2*a+b,2),F(3*b,2)) for a,b in AB]
RAW_Q=[(F(2),F(1)),(F(3,2),F(5,2)),(F(1,2),F(7,2)),(F(-1,2),F(7,2)),(F(-3,2),F(5,2)),(F(-2),F(1)),(F(-2),F(-1)),(F(-3,2),F(-5,2)),(F(-1,2),F(-7,2)),(F(1,2),F(-7,2)),(F(3,2),F(-5,2)),(F(2),F(-1))]
def transform(p):
 x,y=p;return ((6*x+y)/13,(2*y-x)/13)
Q=list(map(transform,RAW_Q));C=list(map(transform,RAW))
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def dot(a,b):return a[0]*b[0]+3*a[1]*b[1]
def norm2(a):return dot(a,a)
def det(a,b):return a[0]*b[1]-a[1]*b[0] # Cartesian determinant divided by sqrt(3).
def orient(a,b,c):return det(sub(b,a),sub(c,a))
FACES=[fa for fa in combinations(range(N),3) if all(norm2(sub(C[i],C[j]))==3*T for i,j in combinations(fa,2))]
ACTIVE=[(1,16,18),(3,13,16),(5,12,13),(7,12,14),(9,14,17),(11,17,18)]
P=[tuple(sum(C[i][d] for i in fa)/3 for d in (0,1)) for fa in FACES]
Z=Q+C+P
EDGES=[(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+FACES.index(fa),B+c) for fa in ACTIVE for c in fa]
WEIGHTS=[F(1,52) if i%2==0 else F(1,39) for i in range(B)]+[F(1,39) if i%2==0 else F(1,52) for i in range(B)]+[F(1,39)]*18
MU=F(-1,156)
ACTIVE_VERTICES=sorted(set(v for e in EDGES for v in e))
