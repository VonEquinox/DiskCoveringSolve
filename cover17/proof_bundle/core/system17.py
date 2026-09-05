"""Explicit 157-variable integer quadratic KKT system for the 17-disk candidate."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=12; N=17; G=82; P=83; TID=82
ACTIVE=[(0, 1, 16), (0, 11, 16), (1, 15, 16), (2, 3, 15), (3, 4, 15), (4, 14, 15), (5, 6, 14), (6, 7, 13), (8, 9, 12), (8, 12, 13), (9, 10, 12), (11, 12, 16), (13, 14, 16)]
FACES=[(0, 1, 16), (0, 11, 16), (1, 2, 15), (1, 15, 16), (2, 3, 15), (3, 4, 15), (4, 5, 14), (4, 14, 15), (5, 6, 14), (6, 7, 13), (6, 13, 14), (7, 8, 13), (8, 9, 12), (8, 12, 13), (9, 10, 12), (10, 11, 12), (11, 12, 16), (12, 13, 16), (13, 14, 16), (14, 15, 16)]
EDGES=[(0, 12), (1, 13), (2, 14), (3, 15), (4, 16), (5, 17), (6, 18), (7, 19), (8, 20), (9, 21), (10, 22), (11, 23), (0, 13), (1, 14), (2, 15), (3, 16), (4, 17), (5, 18), (6, 19), (7, 20), (8, 21), (9, 22), (10, 23), (11, 12), (29, 12), (29, 13), (29, 28), (30, 12), (30, 23), (30, 28), (31, 13), (31, 27), (31, 28), (32, 14), (32, 15), (32, 27), (33, 15), (33, 16), (33, 27), (34, 16), (34, 26), (34, 27), (35, 17), (35, 18), (35, 26), (36, 18), (36, 19), (36, 25), (37, 20), (37, 21), (37, 24), (38, 20), (38, 24), (38, 25), (39, 21), (39, 22), (39, 24), (40, 23), (40, 24), (40, 28), (41, 25), (41, 26), (41, 28)]
M=63; C=74; D=157

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
    d=json.load(open(R/'root17.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
