"""Explicit 148-variable integer quadratic KKT system for the 15-disk candidate."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=11; N=15; G=76; P=77; TID=76
ACTIVE=[(0, 1, 11), (2, 11, 12), (3, 4, 12), (4, 5, 12), (5, 6, 13), (5, 12, 13), (6, 7, 13), (7, 8, 14), (7, 13, 14), (8, 9, 14), (10, 11, 14), (11, 12, 13), (11, 13, 14)]
FACES=[(0, 1, 11), (0, 10, 11), (1, 2, 11), (2, 3, 12), (2, 11, 12), (3, 4, 12), (4, 5, 12), (5, 6, 13), (5, 12, 13), (6, 7, 13), (7, 8, 14), (7, 13, 14), (8, 9, 14), (9, 10, 14), (10, 11, 14), (11, 12, 13), (11, 13, 14)]
EDGES=[(0, 11), (1, 12), (2, 13), (3, 14), (4, 15), (5, 16), (6, 17), (7, 18), (8, 19), (9, 20), (10, 21), (0, 12), (1, 13), (2, 14), (3, 15), (4, 16), (5, 17), (6, 18), (7, 19), (8, 20), (9, 21), (10, 11), (26, 11), (26, 12), (26, 22), (27, 13), (27, 22), (27, 23), (28, 14), (28, 15), (28, 23), (29, 15), (29, 16), (29, 23), (30, 16), (30, 17), (30, 24), (31, 16), (31, 23), (31, 24), (32, 17), (32, 18), (32, 24), (33, 18), (33, 19), (33, 25), (34, 18), (34, 24), (34, 25), (35, 19), (35, 20), (35, 25), (36, 21), (36, 22), (36, 25), (37, 22), (37, 23), (37, 24), (38, 22), (38, 24), (38, 25)]
M=61; C=71; D=148

def build():
    cons=[]
    for u,v in EDGES:
        H={}; b={TID:-1}; c=0
        for d in (0,1):
            a={}; c0=0
            for n,s in [(u,1),(v,-1)]:
                if n: a[2*(n-1)+d]=s
                else: c0+=s*(d==0)
            c+=c0*c0
            for i,x in a.items():
                b[i]=b.get(i,0)+2*c0*x
                for j,y in a.items(): H[i,j]=H.get((i,j),0)+2*x*y
        cons.append((H,b,c))
    for i in range(1,B):
        cons.append(({(2*(i-1)+d,2*(i-1)+d):2 for d in (0,1)},{},-1))
    return cons
CONS=build()

def evaluate(z,zero=0):
    cv=[]; J=[[zero]*D for _ in range(D)]; fv=[zero]*P; fv[TID]=1
    for i,(H,b,c) in enumerate(CONS):
        val=c+sum(v*z[j] for j,v in b.items())
        gr=[b.get(j,zero) for j in range(P)]
        for (j,k),v in H.items():
            val+=v*z[j]*z[k]/2; gr[j]+=v*z[k]
            J[C+j][k]+=z[P+i]*v
        cv.append(val)
        for j,v in enumerate(gr):
            J[i][j]=v; J[C+j][P+i]=v; fv[j]+=z[P+i]*v
    return cv+fv,J

def load():
    d=json.load(open(R/'root15.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
