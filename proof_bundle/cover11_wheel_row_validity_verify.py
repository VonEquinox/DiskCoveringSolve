#!/usr/bin/env python3
"""Pure-Fraction verifier for the wheel circulant-row validity certificate.
No floating or transcendental computation participates in verification.
"""
from fractions import Fraction as F
import json
B=10
C=json.load(open('/mnt/data/cover11_wheel_row_validity.json'))
D=json.load(open('/mnt/data/cover11_wheel_cert_T14438.json'))
QBITS=int(C['qbits']);Q=1<<QBITS
nums=list(map(int,C['mass_nums']));assert C['mass_den']==str(Q)
mass=[F(n,Q) for n in nums];assert len(mass)==3 and min(mass)>0 and sum(mass,F(0))==1
# Graph nodes: anchors 0..9; boundary centers 10..19; center 20; witnesses 21..30.
edges=[];types=[]
for i in range(B):
 edges += [(i,B+i),(i,B+(i+1)%B)];types += [0,0]
for i in range(B):
 p=21+i
 for v in (i,(i+1)%B,10):
  edges.append((p,B+v));types.append(1 if v<10 else 2)
assert len(edges)==50
weights=[mass[0]/20 if t==0 else mass[1]/20 if t==1 else mass[2]/10 for t in types]
assert min(weights)>0 and sum(weights,F(0))==1
adj={i:{} for i in range(31)}
def add(a,b,z):
 if a==b or z==0:return
 adj[a][b]=adj[a].get(b,F(0))+z;adj[b][a]=adj[b].get(a,F(0))+z
for (a,b),w in zip(edges,weights):add(a,b,w)
for v in range(B,31):
 nei=[(u,z) for u,z in adj.get(v,{}).items() if z];s=sum((z for _,z in nei),F(0));assert s>0
 for k,(a,x) in enumerate(nei):
  for b,y in nei[k+1:]:add(a,b,x*y/s)
 for a,_ in nei:adj[a].pop(v,None)
 del adj[v]
conduct=[]
for d in range(1,6):
 vals=[]
 for i in range(B):
  j=(i+d)%B;a,b=(i,j) if i<j else (j,i);vals.append(adj[a].get(b,F(0)))
 assert all(z==vals[0] for z in vals);conduct.append(vals[0])
stored=[F(int(a),int(b)) for a,b in C['conductances']]
assert conduct==stored
A=[F(*D['acoef'][str(d)]) for d in range(1,6)]
assert C['target_acoef']==D['acoef']
# Exact Laplacian difference L(conduct)-L(A).
M=[[F(0) for _ in range(B)] for _ in range(B)]
for i in range(B):
 for j in range(i+1,B):
  d=min(j-i,B-(j-i));z=conduct[d-1]-A[d-1]
  M[i][i]+=z;M[j][j]+=z;M[i][j]-=z;M[j][i]-=z
for row in M:assert sum(row,F(0))==0
# Gauge x_9=0: positive definite 9x9 principal block implies PSD modulo constants.
N=[row[:9] for row in M[:9]];L=[[F(0)]*9 for _ in range(9)];piv=[]
for i in range(9):
 z=N[i][i]-sum((L[i][k]*L[i][k]*piv[k] for k in range(i)),F(0));assert z>0
 piv.append(z);L[i][i]=1
 for j in range(i+1,9):L[j][i]=(N[j][i]-sum((L[j][k]*L[i][k]*piv[k] for k in range(i)),F(0)))/z
stored_piv=[F(int(a),int(b)) for a,b in C['ldlt_pivots']];assert piv==stored_piv
print('PURE-FRACTION WHEEL ROW VALIDITY VERIFIED')
print('minimum exact LDL pivot',float(min(piv)))
