"""Explicit 198-variable integer quadratic KKT system for n=20."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=13; N=20; G=102; P=103; TID=102
ACTIVE=[[0, 15, 16], [1, 2, 16], [3, 16, 17], [4, 5, 17], [5, 17, 18], [6, 7, 18], [7, 8, 19], [7, 18, 19], [8, 9, 19], [9, 10, 14], [9, 14, 19], [11, 12, 15], [11, 14, 15], [13, 14, 15], [13, 14, 19], [13, 15, 16], [13, 16, 17], [13, 17, 18], [13, 18, 19]]
FACES=[[0, 1, 16], [0, 12, 15], [0, 15, 16], [1, 2, 16], [2, 3, 16], [3, 4, 17], [3, 16, 17], [4, 5, 17], [5, 6, 18], [5, 17, 18], [6, 7, 18], [7, 8, 19], [7, 18, 19], [8, 9, 19], [9, 10, 14], [9, 14, 19], [10, 11, 14], [11, 12, 15], [11, 14, 15], [13, 14, 15], [13, 14, 19], [13, 15, 16], [13, 16, 17], [13, 17, 18], [13, 18, 19]]
EDGES=[[0, 13], [1, 14], [2, 15], [3, 16], [4, 17], [5, 18], [6, 19], [7, 20], [8, 21], [9, 22], [10, 23], [11, 24], [12, 25], [0, 14], [1, 15], [2, 16], [3, 17], [4, 18], [5, 19], [6, 20], [7, 21], [8, 22], [9, 23], [10, 24], [11, 25], [12, 13], [33, 13], [33, 28], [33, 29], [34, 14], [34, 15], [34, 29], [35, 16], [35, 29], [35, 30], [36, 17], [36, 18], [36, 30], [37, 18], [37, 30], [37, 31], [38, 19], [38, 20], [38, 31], [39, 20], [39, 21], [39, 32], [40, 20], [40, 31], [40, 32], [41, 21], [41, 22], [41, 32], [42, 22], [42, 23], [42, 27], [43, 22], [43, 27], [43, 32], [44, 24], [44, 25], [44, 28], [45, 24], [45, 27], [45, 28], [46, 26], [46, 27], [46, 28], [47, 26], [47, 27], [47, 32], [48, 26], [48, 28], [48, 29], [49, 26], [49, 29], [49, 30], [50, 26], [50, 30], [50, 31], [51, 26], [51, 31], [51, 32]]
M=83; C=95; D=198

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
    d=json.load(open(R/'root20.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
