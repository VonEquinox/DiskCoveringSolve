"""Canonical 116-variable integer quadratic KKT system for the 14-disk candidate."""
import json
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
B=10; N=14; G=62; P=63; TID=62
ACTIVE=[(1,2,11),(1,10,11),(3,4,11),(4,11,12),(6,7,13),(6,12,13),(8,9,13),(9,10,13)]
FACES=[(0,1,10),(0,9,10),(1,2,11),(1,10,11),(2,3,11),(3,4,11),(4,5,12),(4,11,12),(5,6,12),(6,7,13),(6,12,13),(7,8,13),(8,9,13),(9,10,13),(10,11,12),(10,12,13)]
EDGES=[(i,10+i) for i in range(10)]+[(i,10+(i+1)%10) for i in range(10)]+[(24+j,10+i) for j,f in enumerate(ACTIVE) for i in f]
M=44; C=53; D=116

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
    d=json.load(open(R/'root14.json'));q=int(d['Qx']);yq=int(d['Qy'])
    return d,[F(int(x),q) for x in d['xnum']],[[F(int(x),yq) for x in row] for row in d['Ynum']],F(d['rho'])
