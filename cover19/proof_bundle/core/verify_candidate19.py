"""Exact upper coverage and global candidate-graph lower bound; no root oracle."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
from itertools import combinations
from math import isqrt
import json,time,hashlib
import geometry19 as G
if not __debug__:raise RuntimeError('Assertions must be enabled')
R=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ldl(H):
 n=len(H);L=[[F(i==j) for j in range(n)] for i in range(n)];d=[]
 for i in range(n):
  v=H[i][i]-sum(L[i][k]**2*d[k] for k in range(i));assert v>0,('nonpositive pivot',i,v);d.append(v)
  for j in range(i+1,n):L[j][i]=(H[j][i]-sum(L[j][k]*L[i][k]*d[k] for k in range(i)))/v
 assert all(H[i][j]==sum(L[i][k]*d[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
 return d

def verify():
 outp=R/'candidate19_verified.json';outp.unlink(missing_ok=True);start=time.time()
 assert G.B==12 and G.N==19 and G.T==F(1,13)
 assert set(G.AB)=={(a,b) for a in range(-2,3) for b in range(-2,3) if max(abs(a),abs(b),abs(a+b))<=2}
 assert len(G.FACES)==24 and len(G.ACTIVE)==6 and len(G.EDGES)==42 and len(G.ACTIVE_VERTICES)==36
 assert G.Q[0]==(1,0) and all(G.norm2(q)==1 for q in G.Q)
 assert len(set(G.Q))==12
 # Every nonincident anchor strictly on the left of every polygon edge.
 assert all(G.orient(G.Q[i],G.Q[(i+1)%12],G.Q[j])>0 for i in range(12) for j in range(12) if j not in (i,(i+1)%12))
 # Each center triangle is equilateral and has the stated common point.
 for fa,p in zip(G.FACES,G.P):
  assert all(G.norm2(G.sub(p,G.C[i]))==G.T for i in fa)
 # Oriented covering chain, including six explicitly allowed degenerate faces.
 tris=[]
 for fa in G.FACES:
  t=tuple(12+i for i in fa)
  if G.orient(*(G.Z[i] for i in t))<0:t=(t[0],t[2],t[1])
  assert G.orient(*(G.Z[i] for i in t))>0;tris.append(t)
 for i in range(12):
  a,b,c=G.Q[(i-1)%12],G.Q[i],G.C[i]
  assert G.norm2(G.sub(a,c))==G.T==G.norm2(G.sub(b,c))
  assert G.det(a,b)>0 and G.dot(a,b)>0 and G.det(a,c)>0 and G.det(c,b)>0
  # First collar triangle has common point q_i; second lies in disk c_i.
  assert G.norm2(G.sub(G.Q[i],G.C[(i+1)%12]))==G.T
  tris.extend([(12+i,i,12+(i+1)%12),((i-1)%12,i,12+i)])
 assert len(tris)==48
 ori=[G.orient(*(G.Z[i] for i in t)) for t in tris];assert min(ori)==0 and sum(v==0 for v in ori)==6
 chain=Counter()
 for a,b,c in tris:
  for u,v in [(a,b),(b,c),(c,a)]:chain[u,v]+=1;chain[v,u]-=1
 expected=Counter()
 for i in range(12):expected[i,(i+1)%12]+=1;expected[(i+1)%12,i]-=1
 assert {k:v for k,v in chain.items() if v}=={k:v for k,v in expected.items() if v}
 # Exact equilibrium, normalization, and equal length on the positive subgraph.
 assert len(G.WEIGHTS)==len(G.EDGES) and all(w>0 for w in G.WEIGHTS) and sum(G.WEIGHTS)==1
 forces=[[F(0),F(0)] for _ in G.Z]
 for (u,v),w in zip(G.EDGES,G.WEIGHTS):
  d=G.sub(G.Z[u],G.Z[v]);assert G.norm2(d)==G.T
  for k in range(2):forces[u][k]+=w*d[k];forces[v][k]-=w*d[k]
 for i in range(12):
  for k in range(2):assert forces[i][k]+G.MU*G.Q[i][k]==0
 assert all(forces[i]==[0,0] for i in G.ACTIVE_VERTICES if i>=12)
 # Exact star-mesh elimination of all free positive-graph nodes.
 C={v:{} for v in G.ACTIVE_VERTICES}
 for (u,v),w in zip(G.EDGES,G.WEIGHTS):C[u][v]=C[u].get(v,F(0))+w;C[v][u]=C[v].get(u,F(0))+w
 elimination=[]
 for v in [i for i in G.ACTIVE_VERTICES if i>=12]:
  nei=list(C[v]);s=sum(C[v].values());assert s>0;elimination.append(s)
  for u,w in combinations(nei,2):
   d=C[v][u]*C[v][w]/s;C[u][w]=C[u].get(w,F(0))+d;C[w][u]=C[w].get(u,F(0))+d
  for u in nei:del C[u][v]
  del C[v]
 assert set(C)==set(range(12))
 conduct=[(i,j,C[i].get(j,F(0))) for i in range(12) for j in range(i+1,12)]
 assert all(c>=0 for i,j,c in conduct)
 # Independent closed form from a six-node cycle with scaled matrix
 # 10 I - 3(S+S^T), whose eigenvalues are 4,7,13,16,13,7.
 cos6=[F(1),F(1,2),F(-1,2),F(-1),F(-1,2),F(1,2)]
 gg=[(F(1,4)+F(2,7)*cos6[d%6]+F(2,13)*cos6[(2*d)%6]+F((-1)**d,16))/351 for d in range(4)]
 for i,j,c in conduct:
  a,b=i//2,j//2;dd=min((b-a)%6,(a-b)%6)
  expected=gg[dd]+(F(1,117) if a==b else 0)
  if (j==i+1 and j%2==0) or (i==0 and j==11):expected+=F(1,104)
  assert c==expected
 assert sum(c*G.norm2(G.sub(G.Q[i],G.Q[j])) for i,j,c in conduct)==G.T
 # delta=2 asin(1/sqrt13), cos(delta)=11/13. For k<=5, k delta<pi.
 # cos(pi/5)=(1+sqrt5)/4 <11/13, because sqrt5<31/13.
 assert F(31,13)**2>5
 cheb=[F(1),F(11,13)]
 for k in range(2,6):cheb.append(2*F(11,13)*cheb[-1]-cheb[-2])
 cheb.append(F(-1))
 H=[[F(0) for _ in range(11)] for _ in range(11)]
 for i,j,c in conduct:
  k=min(j-i,12-j+i);v=2*c*cheb[k]
  if i:H[i-1][i-1]+=v
  H[j-1][j-1]+=v
  if i:H[i-1][j-1]-=v;H[j-1][i-1]-=v
 piv=ldl(H)
 # A smaller independent certificate: a circulant lower matrix on all 12 angles.
 weights=[None]+[min(2*c*cheb[k] for i,j,c in conduct if min(j-i,12-j+i)==k) for k in range(1,7)]
 cosine=[(F(1),F(0)),(F(0),F(1,2)),(F(1,2),F(0)),(F(0),F(0)),(F(-1,2),F(0)),(F(0),F(-1,2)),(F(-1),F(0)),(F(0),F(-1,2)),(F(-1,2),F(0)),(F(0),F(0)),(F(1,2),F(0)),(F(0),F(1,2))]
 sqrtlo,sqrthi=F(173205,100000),F(173206,100000)
 assert sqrtlo**2<3<sqrthi**2
 spectrum=[]
 for m in range(1,12):
  rational=F(0);radical=F(0)
  for k in range(1,7):
   ca,cb=cosine[(m*k)%12];v=(2 if k<6 else 1)*weights[k]
   rational+=v*(1-ca);radical-=v*cb
  lower=rational+radical*(sqrtlo if radical>=0 else sqrthi)
  assert lower>F(7,5000)
  spectrum.append({'mode':m,'rational':str(rational),'sqrt3_coefficient':str(radical),'lower':str(lower)})
 # Independently verify angular stationarity after harmonic elimination.
 for i in range(12):assert sum(C[i].get(j,F(0))*G.det(G.Q[i],G.Q[j]) for j in range(12))==0
 digits=60;Q=10**digits;rlo=isqrt(Q*Q//13);Rlo=isqrt(13*Q*Q)
 assert 13*rlo*rlo<Q*Q<13*(rlo+1)**2 and Rlo*Rlo<13*Q*Q<(Rlo+1)**2
 dec=lambda n:str(n//Q)+'.'+str(n%Q).zfill(digits)
 out={'verified':True,'exact_radius_squared':'1/13','r_interval':[dec(rlo),dec(rlo+1)],'R_interval':[dec(Rlo),dec(Rlo+1)],'full_centers':19,'boundary_anchors':12,'center_faces':24,'covering_chain_faces':48,'positive_triangles':42,'degenerate_triangles':6,'boundary_caps':12,'positive_graph_vertices':36,'positive_rods':42,'radial_multiplier':str(G.MU),'effective_conductances':[[i,j,str(c)] for i,j,c in conduct],'hessian_lower':[[str(v) for v in row] for row in H],'ldl_pivots':list(map(str,piv)),'minimum_ldl_pivot':str(min(piv)),'global_strict_convexity_verified':True,'circulant_lower_weights':list(map(str,weights[1:])),'fourier_eigenvalues':spectrum,'nonconstant_spectral_gap_lower':'7/5000','geometry_sha256':sha(R/'geometry19.py'),'seconds':time.time()-start}
 outp.write_text(json.dumps(out,indent=2));print('EXACT CANDIDATE UPPER + GLOBAL LOWER VERIFIED',out['r_interval'],'min pivot',float(min(piv)),flush=True);return out
if __name__=='__main__':verify()
