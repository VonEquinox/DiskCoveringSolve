#!/usr/bin/env python3
"""Exact validity certificate for the circulant wheel energy row.

The old wheel branch tree stores only its five circulant anchor
conductances.  This script supplies a self-contained graph witness: a positive
300-bit rational symmetric edge weighting of the 50-edge wheel graph.  Its
exact Dirichlet-to-Neumann Laplacian dominates the stored circulant Laplacian
in Loewner order.  Therefore every stored wheel leaf uses a genuine lower
bound for the graph minimax radius.
"""
from fractions import Fraction as F
import json,mpmath as mp
mp.mp.dps=200
B=10
# Graph nodes: anchors 0..9, centers 10..20, face witnesses 21..30.
edges=[];types=[]
for i in range(B):
    edges += [(i,B+i),(i,B+(i+1)%B)];types += [0,0]
for i in range(B):
    p=21+i
    for v in (i,(i+1)%B,10):
        edges.append((p,B+v));types.append(1 if v<10 else 2)
assert len(edges)==50
# Rationalized golden-ratio symmetric type masses.  The formula is only a
# proposal; all proof checks below use the resulting exact rational integers.
QBITS=300;Q=1<<QBITS;s=mp.sqrt(5)
real=[mp.mpf(2)/(s+2),(s-1)/(s+2),mp.mpf(1)/(s+2)]
nums=[int(mp.floor(z*Q)) for z in real];nums[2]+=Q-sum(nums)
mass=[F(n,Q) for n in nums]
assert min(mass)>0 and sum(mass,F(0))==1
weights=[mass[0]/20 if t==0 else mass[1]/20 if t==1 else mass[2]/10 for t in types]
assert min(weights)>0 and sum(weights,F(0))==1
# Exact star-mesh elimination.
nnode=31;adj={i:{} for i in range(nnode)}
def add(a,b,z):
    if a==b or z==0:return
    adj[a][b]=adj[a].get(b,F(0))+z;adj[b][a]=adj[b].get(a,F(0))+z
for (a,b),w in zip(edges,weights):add(a,b,w)
for v in range(B,nnode):
    nei=[(u,z) for u,z in adj.get(v,{}).items() if z]
    ss=sum((z for _,z in nei),F(0));assert ss>0
    for k,(a,x) in enumerate(nei):
        for b,y in nei[k+1:]:add(a,b,x*y/ss)
    for a,_ in nei:adj[a].pop(v,None)
    del adj[v]
# Symmetry is checked exactly, not assumed.
C=[]
for d in range(1,6):
    vals=[]
    for i in range(B):
        j=(i+d)%B;a,b=(i,j) if i<j else (j,i);vals.append(adj[a].get(b,F(0)))
    assert all(z==vals[0] for z in vals)
    C.append(vals[0])
D=json.load(open('/mnt/data/cover11_wheel_cert_T14438.json'))
A=[F(*D['acoef'][str(d)]) for d in range(1,6)]
# Laplacian difference L(C)-L(A).
M=[[F(0) for _ in range(B)] for _ in range(B)]
for i in range(B):
    for j in range(i+1,B):
        d=min(j-i,B-(j-i));z=C[d-1]-A[d-1]
        M[i][i]+=z;M[j][j]+=z;M[i][j]-=z;M[j][i]-=z
for i in range(B):assert sum(M[i],F(0))==0
# Because the constant vector is the kernel, set x_9=0.  Exact LDL^T of the
# leading 9x9 block proves positive definiteness on the quotient by constants.
N=[row[:9] for row in M[:9]];L=[[F(0)]*9 for _ in range(9)];piv=[]
for i in range(9):
    z=N[i][i]-sum((L[i][k]*L[i][k]*piv[k] for k in range(i)),F(0))
    assert z>0,(i,z);piv.append(z);L[i][i]=1
    for j in range(i+1,9):
        L[j][i]=(N[j][i]-sum((L[j][k]*L[i][k]*piv[k] for k in range(i)),F(0)))/z
out={
 'qbits':QBITS,'mass_nums':[str(n) for n in nums],'mass_den':str(Q),
 'conductances':[[str(z.numerator),str(z.denominator)] for z in C],
 'target_acoef':D['acoef'],
 'ldlt_pivots':[[str(z.numerator),str(z.denominator)] for z in piv],
 'min_pivot_float':float(min(piv)),
}
open('/mnt/data/cover11_wheel_row_validity.json','w').write(json.dumps(out,separators=(',',':')))
print('WHEEL CIRCULANT ROW VALIDITY VERIFIED')
print('type masses',[float(z) for z in mass])
print('coefficient differences',[float(c-a) for c,a in zip(C,A)])
print('minimum exact LDL pivot',float(min(piv)))
