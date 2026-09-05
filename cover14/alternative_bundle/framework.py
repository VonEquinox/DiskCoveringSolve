"""Exact canonical 14-disk candidate framework. All polynomial coefficients are integers."""
from pathlib import Path
N_DISKS=14;B=10
ALL_FACES=[(0,1,10),(0,9,10),(1,2,10),(2,3,11),(2,10,11),(3,4,11),(4,5,12),(4,11,12),(5,6,12),(6,7,12),(7,8,13),(7,12,13),(8,9,13),(9,10,13),(10,11,13),(11,12,13)]
ACTIVE_IDS=[1,2,4,6,7,9,11,13]
ACTIVE_FACES=[ALL_FACES[k] for k in ACTIVE_IDS]
NN=B+N_DISKS+len(ACTIVE_FACES)
ND=2*(NN-1);TID=ND
EDGES=[(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N_DISKS+j,B+c) for j,f in enumerate(ACTIVE_FACES) for c in f]
NE=len(EDGES);NC=NE+B-1;WSTART=ND+1;MUSTART=WSTART+NE;N=ND+1+NC
assert (ND,NE,N)==(62,44,116)

def coord(node,d):return (2*(node-1)+d) if node else None

def build():
 # Each constraint is (sparse symmetric Hessian H, linear b, constant c).
 cons=[]
 for u,v in EDGES:
  H={};b={TID:-1};cc=0
  for d in range(2):
   a={};c=0
   for node,sgn in ((u,1),(v,-1)):
    if node:a[coord(node,d)]=sgn
    else:c+=sgn*(1 if d==0 else 0)
   cc+=c*c
   for i,aa in a.items():
    b[i]=b.get(i,0)+2*c*aa
    for j,bb in a.items():H[i,j]=2*aa*bb
  cons.append((H,{k:v for k,v in b.items() if v},cc))
 for node in range(1,B):cons.append(({(coord(node,d),coord(node,d)):2 for d in range(2)},{},-1))
 return cons
CONS=build()

def values_jac(z):
 """Works for Python floats, fractions.Fraction and mpmath.mpf."""
 zero=z[0]*0;one=zero+1
 x=z[:ND+1];lam=z[WSTART:]
 G=[];vals=[]
 for H,b,c in CONS:
  val=zero+c
  for j,v in b.items():val+=v*x[j]
  # Integral Hessians have even entries: division by 2 is exact in Z.
  for (i,j),v in H.items():val+=(v//2)*x[i]*x[j]
  vals.append(val)
  g=[zero]*(ND+1)
  for j,v in b.items():g[j]+=v
  for (i,j),v in H.items():g[i]+=v*x[j]
  G.append(g)
 stat=[(one if j==TID else zero)+sum((lam[k]*G[k][j] for k in range(NC)),zero) for j in range(ND+1)]
 J=[[zero]*N for _ in range(N)]
 for i in range(NC):
  for j in range(ND+1):J[i][j]=G[i][j];J[NC+j][WSTART+i]=G[i][j]
  for (j,k),v in CONS[i][0].items():J[NC+j][k]+=lam[i]*v
 return vals+stat,J

def j_lipschitz_rows():
 # For common root box radius rho, ||J_row(z)-J_row(z0)||_1 <= rho*D[row].
 D=[0]*N
 for i,(H,b,c) in enumerate(CONS):
  D[i]=sum(abs(v) for v in H.values())
  for (j,k),v in H.items():D[NC+j]+=2*abs(v)
 return D
