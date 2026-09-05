"""Explicit 174-variable integer quadratic KKT system for n=18."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=12; N=18; G=90; P=91; TID=90
ACTIVE=[[0, 1, 14], [1, 2, 15], [1, 14, 15], [2, 3, 15], [3, 4, 16], [3, 15, 16], [5, 16, 17], [6, 7, 17], [7, 8, 17], [8, 9, 13], [8, 12, 13], [8, 12, 17], [9, 10, 13], [11, 13, 14], [12, 14, 15], [12, 15, 16]]
FACES=[[0, 1, 14], [0, 11, 14], [1, 2, 15], [1, 14, 15], [2, 3, 15], [3, 4, 16], [3, 15, 16], [4, 5, 16], [5, 6, 17], [5, 16, 17], [6, 7, 17], [7, 8, 17], [8, 9, 13], [8, 12, 13], [8, 12, 17], [9, 10, 13], [10, 11, 13], [11, 13, 14], [12, 13, 14], [12, 14, 15], [12, 15, 16], [12, 16, 17]]
EDGES=[[0, 12], [1, 13], [2, 14], [3, 15], [4, 16], [5, 17], [6, 18], [7, 19], [8, 20], [9, 21], [10, 22], [11, 23], [0, 13], [1, 14], [2, 15], [3, 16], [4, 17], [5, 18], [6, 19], [7, 20], [8, 21], [9, 22], [10, 23], [11, 12], [30, 12], [30, 13], [30, 26], [31, 13], [31, 14], [31, 27], [32, 13], [32, 26], [32, 27], [33, 14], [33, 15], [33, 27], [34, 15], [34, 16], [34, 28], [35, 15], [35, 27], [35, 28], [36, 17], [36, 28], [36, 29], [37, 18], [37, 19], [37, 29], [38, 19], [38, 20], [38, 29], [39, 20], [39, 21], [39, 25], [40, 20], [40, 24], [40, 25], [41, 20], [41, 24], [41, 29], [42, 21], [42, 22], [42, 25], [43, 23], [43, 25], [43, 26], [44, 24], [44, 26], [44, 27], [45, 24], [45, 27], [45, 28]]
M=72; C=83; D=174

ACTIVE=[tuple(f) for f in ACTIVE]
FACES=[tuple(f) for f in FACES]
EDGES=[tuple(e) for e in EDGES]

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
