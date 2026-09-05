"""Explicit 174-variable integer quadratic KKT system for 18 disks."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=12; N=18; G=90; P=91; TID=90
ACTIVE=[(0, 1, 14), (0, 11, 14), (1, 2, 13), (1, 13, 14), (3, 12, 13), (4, 5, 12), (5, 6, 12), (6, 7, 16), (6, 12, 17), (6, 16, 17), (7, 8, 16), (9, 15, 16), (10, 11, 15), (11, 14, 15), (13, 14, 17), (14, 15, 17)]
FACES=[(0, 1, 14), (0, 11, 14), (1, 2, 13), (1, 13, 14), (2, 3, 13), (3, 4, 12), (3, 12, 13), (4, 5, 12), (5, 6, 12), (6, 7, 16), (6, 12, 17), (6, 16, 17), (7, 8, 16), (8, 9, 16), (9, 10, 15), (9, 15, 16), (10, 11, 15), (11, 14, 15), (12, 13, 17), (13, 14, 17), (14, 15, 17), (15, 16, 17)]
EDGES=[(0, 12), (1, 13), (2, 14), (3, 15), (4, 16), (5, 17), (6, 18), (7, 19), (8, 20), (9, 21), (10, 22), (11, 23), (0, 13), (1, 14), (2, 15), (3, 16), (4, 17), (5, 18), (6, 19), (7, 20), (8, 21), (9, 22), (10, 23), (11, 12), (30, 12), (30, 13), (30, 26), (31, 12), (31, 23), (31, 26), (32, 13), (32, 14), (32, 25), (33, 13), (33, 25), (33, 26), (34, 15), (34, 24), (34, 25), (35, 16), (35, 17), (35, 24), (36, 17), (36, 18), (36, 24), (37, 18), (37, 19), (37, 28), (38, 18), (38, 24), (38, 29), (39, 18), (39, 28), (39, 29), (40, 19), (40, 20), (40, 28), (41, 21), (41, 27), (41, 28), (42, 22), (42, 23), (42, 27), (43, 23), (43, 26), (43, 27), (44, 25), (44, 26), (44, 29), (45, 26), (45, 27), (45, 29)]
M=72; C=83; D=174

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
    d=json.load(open(R/'root18.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
