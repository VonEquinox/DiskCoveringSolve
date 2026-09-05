"""Explicit 140-variable integer quadratic KKT system for the 16-disk candidate."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=12; N=16; G=74; P=75; TID=74
ACTIVE=[(0, 1, 12), (2, 12, 13), (3, 4, 13), (5, 13, 14), (6, 7, 14), (8, 14, 15), (9, 10, 15), (11, 12, 15), (12, 13, 15), (13, 14, 15)]
FACES=[(0, 1, 12), (0, 11, 12), (1, 2, 12), (2, 3, 13), (2, 12, 13), (3, 4, 13), (4, 5, 13), (5, 6, 14), (5, 13, 14), (6, 7, 14), (7, 8, 14), (8, 9, 15), (8, 14, 15), (9, 10, 15), (10, 11, 15), (11, 12, 15), (12, 13, 15), (13, 14, 15)]
EDGES=[(0, 12), (1, 13), (2, 14), (3, 15), (4, 16), (5, 17), (6, 18), (7, 19), (8, 20), (9, 21), (10, 22), (11, 23), (0, 13), (1, 14), (2, 15), (3, 16), (4, 17), (5, 18), (6, 19), (7, 20), (8, 21), (9, 22), (10, 23), (11, 12), (28, 12), (28, 13), (28, 24), (29, 14), (29, 24), (29, 25), (30, 15), (30, 16), (30, 25), (31, 17), (31, 25), (31, 26), (32, 18), (32, 19), (32, 26), (33, 20), (33, 26), (33, 27), (34, 21), (34, 22), (34, 27), (35, 23), (35, 24), (35, 27), (36, 24), (36, 25), (36, 27), (37, 25), (37, 26), (37, 27)]
M=54; C=65; D=140

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
    d=json.load(open(R/'root16.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
