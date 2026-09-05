#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
import json,time,math
R=Path(__file__).resolve().parent; CORE=R

def load():
 D=json.load(open(CORE/'cover13_krawczyk_cert.json'))
 L=json.load(open(R/'cover13_core_layout.json'))
 Qx=int(D['Qx']); Qy=int(D['Qy']); rho=F(int(D['rho_num']),int(D['rho_den']))
 z=[F(int(a),Qx) for a in D['xnum']]
 Yold=[[F(int(a),Qy) for a in row] for row in D['Ynum']]
 # Certificate row order: rods; free balances; anchor balances; circle norms;
 # sum(weights)-1. The verifier uses constraints followed by full stationarity.
 # Jold = P Jnew; transform inverse by Ynew = Yold P, EXACTLY.
 permutation=list(range(47))+list(range(74,118))+list(range(56,74))+list(range(47,56))+[118]
 scales=[F(1)]*47+[F(1,2)]*62+[F(1)]*9+[F(-1)]
 assert sorted(permutation)==list(range(119)) and len(scales)==119
 Y=[[F(0) for _ in range(119)] for _ in range(119)]
 for i in range(119):
  for old,new in enumerate(permutation):Y[i][new]=Yold[i][old]*scales[old]
 assert int(D['N'])==119 and len(z)==119 and len(Y)==119 and all(len(row)==119 for row in Y)
 return D,L,z,Y,rho

def idxmap(layout):
 mp={};k=0
 if layout=='qcp':
  for i in range(1,10):mp[('q',i)]=(k,k+1);k+=2
  for i in range(13):mp[('c',i)]=(k,k+1);k+=2
  for i in range(9):mp[('p',i)]=(k,k+1);k+=2
 elif layout=='cqp':
  for i in range(13):mp[('c',i)]=(k,k+1);k+=2
  for i in range(1,10):mp[('q',i)]=(k,k+1);k+=2
  for i in range(9):mp[('p',i)]=(k,k+1);k+=2
 elif layout=='qpc':
  for i in range(1,10):mp[('q',i)]=(k,k+1);k+=2
  for i in range(9):mp[('p',i)]=(k,k+1);k+=2
  for i in range(13):mp[('c',i)]=(k,k+1);k+=2
 else:raise ValueError(layout)
 assert k==62
 return mp

def build(L):
 layout=L['layout']; mp=idxmap(layout); order=list(map(int,L['boundary_order'])); triples=[tuple(map(int,t)) for t in L['active_triples']]
 def affine(pt,d):
  if pt==('q',0):return {},F(1 if d==0 else 0)
  j=mp[pt][d];return {j:F(1)},F(0)
 ee=[]
 for i in range(10):
  if L['boundary_convention']=='prev': ee += [(('q',i),('c',order[(i-1)%10])),(('q',i),('c',order[i]))]
  else:ee += [(('q',i),('c',order[i])),(('q',i),('c',order[(i+1)%10]))]
 tris=list(reversed(triples)) if L['reverse_triples'] else triples
 for pi,t in enumerate(tris):
  for v in t:ee.append((('p',pi),('c',v)))
 assert len(ee)==47
 rods=[]
 for u,v in ee:
  H={};b={62:F(-1)};c=F(0)
  for d in (0,1):
   du,cu=affine(u,d);dv,cv=affine(v,d);a=dict(du)
   for j,x in dv.items():a[j]=a.get(j,F(0))-x
   c0=cu-cv
   for i,x in a.items():
    b[i]=b.get(i,F(0))+2*c0*x
    for j,y in a.items():H[i,j]=H.get((i,j),F(0))+2*x*y
   c+=c0*c0
  rods.append((H,b,c))
 norms=[]
 for qi in range(1,10):
  H={};b={};c=F(-1)
  for d in (0,1):
   a,c0=affine(('q',qi),d)
   for i,x in a.items():
    b[i]=b.get(i,F(0))+2*c0*x
    for j,y in a.items():H[i,j]=H.get((i,j),F(0))+2*x*y
   c+=c0*c0
  norms.append((H,b,c))
 cons=rods+norms if L['constraint_order']=='rn' else norms+rods
 return cons,ee

def val(c,x):
 H,b,cc=c
 return cc+sum(v*x[i] for i,v in b.items())+sum(F(1,2)*a*x[i]*x[j] for (i,j),a in H.items())
def grad(c,x):
 H,b,cc=c;g=[b.get(i,F(0)) for i in range(63)]
 for (i,j),a in H.items():g[i]+=a*x[j]
 return g

def matmul(A,B):
 n=len(A);p=len(B[0]);cols=[[(k,B[k][j]) for k in range(len(B)) if B[k][j]] for j in range(p)]
 return [[sum((A[i][k]*v for k,v in cols[j] if A[i][k]),F(0)) for j in range(p)] for i in range(n)]

def verify():
 t0=time.time();D,L,z,Y,rho=load();cons,edges=build(L);x=z[:63];lam=z[63:];ls=F(int(L['lambda_sign']));os=F(int(L['objective_sign']))
 assert len(cons)==56
 cv=[val(c,x) for c in cons];G=[grad(c,x) for c in cons]
 stat=[(os if j==62 else F(0))+ls*sum(lam[i]*G[i][j] for i in range(56)) for j in range(63)]
 F0=cv+stat
 n=119;J=[[F(0) for _ in range(n)] for _ in range(n)]
 for i in range(56):
  for j in range(63):J[i][j]=G[i][j];J[56+j][63+i]=ls*G[i][j]
 for j in range(63):
  for k in range(63):J[56+j][k]=ls*sum(lam[i]*cons[i][0].get((j,k),F(0)) for i in range(56))
 YF=[sum(Y[i][k]*F0[k] for k in range(n)) for i in range(n)]
 YJ=matmul(Y,J)
 E=[[F(1 if i==j else 0)-YJ[i][j] for j in range(n)] for i in range(n)]
 # coefficient matrix Dcoef with |J(X)-J0| <= rho Dcoef entrywise
 Dc=[[F(0) for _ in range(n)] for _ in range(n)]
 for i,c in enumerate(cons):
  H=c[0]
  for j in range(63):Dc[i][j]=sum(abs(H.get((j,k),F(0))) for k in range(63))
 for j in range(63):
  for k in range(63):Dc[56+j][k]=sum(abs(cons[i][0].get((j,k),F(0))) for i in range(56))
  for i in range(56):Dc[56+j][63+i]=sum(abs(cons[i][0].get((j,k),F(0))) for k in range(63))
 Dr=[sum(row,F(0)) for row in Dc]
 ratios=[];contr=[]
 for i in range(n):
  r=abs(YF[i])/rho+sum(abs(a) for a in E[i])+rho*sum(abs(Y[i][k])*Dr[k] for k in range(n))
  ratios.append(r)
  contr.append(sum(abs(a) for a in E[i])+rho*sum(abs(Y[i][k])*Dr[k] for k in range(n)))
 ratio=max(ratios);q=max(contr)
 assert ratio<1,(float(ratio),ratios.index(ratio));assert q<1
 # Positive physical rod multipliers and normalization from t stationarity.
 physical=[ls*lam[i] for i,l in enumerate(L['constraint_labels']) if l[0]=='rod']
 assert len(physical)==47 and min(v-rho for v in physical)>0
 tsum=F(0)
 for i,l in enumerate(L['constraint_labels']):
  if l[0]=='rod':tsum+=physical[i]
 # derivative g_e/dt=-1; stationarity os - sum physical weights=0
 assert abs(tsum-os)<F(1,10**40)
 ti=62;tlo=x[ti]-rho;thi=x[ti]+rho;assert tlo>0
 out={'verified':True,'inverse_interface_repaired':True,'dimension':n,'box_radius':str(rho),'max_krawczyk_ratio':str(ratio),'max_krawczyk_ratio_float':float(ratio),'contraction':str(q),'contraction_float':float(q),'midpoint_residual_inf':str(max(abs(a) for a in F0)),'min_rod_weight':str(min(v-rho for v in physical)),'min_rod_weight_float':float(min(v-rho for v in physical)),'rod_weight_sum_midpoint':str(tsum),'t_interval':[str(tlo),str(thi)],'r_interval_float':[math.sqrt(float(tlo)),math.sqrt(float(thi))],'layout':{k:L[k] for k in ('layout','boundary_convention','reverse_triples','constraint_order','objective_sign','lambda_sign')},'seconds':time.time()-t0}
 (R/'cover13_krawczyk_core_verified.json').write_text(json.dumps(out,indent=2));(R/'KRAWCZYK_CORE_OK').write_text('ok\n');print(json.dumps(out,indent=2))
if __name__=='__main__':verify()
