#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import json,math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
ROOT=Path(__file__).resolve().parent
D=json.load(open(ROOT/'cover12_symmetric_kkt_cert.json'));Q=int(D['Qx'])
x={n:float(F(int(a),Q)) for n,a in zip(D['names'],D['xnum'])}
r,u,v,w,c,s=[x[k] for k in ['r','u','v','w','c','s']]
z=math.sqrt(3)
def rot(p):
 a,b=p;return np.array((-(a+z*b)/2,(z*a-b)/2))
def orbit(p):
 p=np.array(p,float);return [p,rot(p),rot(rot(p))]
raw=orbit((u*c,u*s))+orbit((u*c,-u*s))+orbit((-v,0))+orbit((w,0))
order=[6,5,4,7,1,0,8,3,2,9,10,11]
C=np.array([raw[i] for i in order])
faces=[(0,1,11),(0,8,10),(0,10,11),(1,2,11),(2,3,11),(3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),(6,9,10),(7,8,10),(9,10,11)]
fig,ax=plt.subplots(figsize=(8,8))
ax.add_patch(Circle((0,0),1,fill=False,linewidth=2))
for p in C:ax.add_patch(Circle(tuple(p),r,fill=False,linewidth=1))
for f in faces:
 p=np.vstack([C[list(f)],C[f[0]]]);ax.plot(p[:,0],p[:,1],linewidth=.6)
ax.scatter(C[:9,0],C[:9,1],s=18)
ax.scatter(C[9:,0],C[9:,1],s=28)
for i,p in enumerate(C):ax.text(p[0]+.015,p[1]+.015,str(i),fontsize=8)
ax.set_aspect('equal');ax.set_xlim(-1.06,1.06);ax.set_ylim(-1.06,1.06)
ax.set_title(r'$n=12$ candidate: $r=0.3611029637445086\ldots$')
ax.set_xlabel('$x$');ax.set_ylabel('$y$')
fig.tight_layout();fig.savefig(ROOT/'cover12_candidate.png',dpi=220);plt.close(fig)
print(ROOT/'cover12_candidate.png')
