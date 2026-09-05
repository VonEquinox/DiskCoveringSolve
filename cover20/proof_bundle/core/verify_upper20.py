"""Exact upper cover for n=20: 51 triangles and 13 circular caps."""
from fractions import Fraction as F
from collections import Counter,defaultdict
from pathlib import Path
import json,time,hashlib
import system20 as S
import interval20 as I
from exact_arcs import sincos_turn,SCALE
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def unit(u):
 co,si=sincos_turn(u)
 return [(F(co[0],SCALE),F(co[1],SCALE)),(F(si[0],SCALE),F(si[1],SCALE))]

def verify():
 st=time.time();d,z,iv,p=I.rootbox();tlo,thi=iv[S.TID];B=S.B;N=S.N
 edges={frozenset(e) for e in S.EDGES};active_w={tuple(f):B+N+i for i,f in enumerate(S.ACTIVE)};slack=[];equality_faces=0
 assert len(S.FACES)==2*N-B-2 and len(set(map(tuple,S.FACES)))==2*N-B-2
 for face in S.FACES:
  if tuple(face) in active_w:
   assert all(frozenset((B+c,active_w[tuple(face)])) in edges for c in face);equality_faces+=1
  else:
   a,b,c=[[z[2*(B+i-1)+j] for j in (0,1)] for i in face]
   u=[2*(b[j]-a[j]) for j in (0,1)];v=[2*(c[j]-a[j]) for j in (0,1)]
   bu=sum(x*x for x in b)-sum(x*x for x in a);bv=sum(x*x for x in c)-sum(x*x for x in a);de=u[0]*v[1]-u[1]*v[0];assert de
   w=[(bu*v[1]-u[1]*bv)/de,(u[0]*bv-bu*v[0])/de];w=[F(round(q*10**12),10**12) for q in w];wiv=[I.point(q) for q in w]
   mar=min(tlo-I.norm2(I.vsub(p[B+i],wiv))[1] for i in face);assert mar>0
   slack.append({'face':list(face),'witness':list(map(str,w)),'margin':str(mar)})
 assert equality_faces==len(S.ACTIVE) and len(slack)==len(S.FACES)-len(S.ACTIVE)
 aa=json.load(open(R/'anchor_angles20.json'));cent=list(map(F,aa['centers']));h=F(aa['halfwidth'])
 assert len(cent)==B-1 and 0<cent[0]-h and cent[-1]+h<1 and h>0
 assert all(cent[i]+h<cent[i+1]-h for i in range(B-2))
 for i,m in enumerate(cent,1):
  l,u=m-h,m+h;assert u-l<F(1,4)
  assert I.cross(unit(l),p[i])[0]>0 and I.cross(p[i],unit(u))[0]>0 and I.dot(unit(m),p[i])[0]>0
 # Every consecutive true arc is strictly less than a quarter turn.
 aa_lo=[F(0)]+[c-h for c in cent]+[F(1)]
 aa_hi=[F(0)]+[c+h for c in cent]+[F(1)]
 assert all(0<aa_lo[i+1]-aa_hi[i] and aa_hi[i+1]-aa_lo[i]<F(1,4) for i in range(B))
 minarc=None
 for i in range(B):
  prev=(i-1)%B;ctr=p[B+i]
  a=I.cross(p[prev],ctr)[0];b=I.cross(ctr,p[i])[0];assert a>0 and b>0 and I.dot(p[prev],p[i])[0]>0
  assert frozenset((prev,B+i)) in edges and frozenset((i,B+i)) in edges
  mar=min(a,b);minarc=mar if minarc is None else min(minarc,mar)
 triangles=[]
 for i in range(B):triangles.extend([((i-1)%B,i,B+i),(i,B+(i+1)%B,B+i)])
 triangles.extend(tuple(B+c for c in face) for face in S.FACES)
 directed=Counter();ec=Counter();links=defaultdict(list);minarea=None;oriented=[]
 for tri in triangles:
  a,b,c=tri;de=I.cross(I.vsub(p[b],p[a]),I.vsub(p[c],p[a]))
  if de[1]<0:b,c=c,b;de=I.cross(I.vsub(p[b],p[a]),I.vsub(p[c],p[a]))
  assert de[0]>0,(tri,de);minarea=de[0] if minarea is None else min(minarea,de[0])
  tri=(a,b,c);oriented.append(list(tri))
  for j in range(3):
   u,v=tri[j],tri[(j+1)%3];directed[u,v]+=1;ec[tuple(sorted((u,v)))]+=1;links[u].append((tri[(j+1)%3],tri[(j+2)%3]))
 assert set(links)==set(range(B+N)) and len(triangles)==2*N+B-2 and len(ec)==3*N+2*B-3 and B+N-len(ec)+len(triangles)==1
 boundary={tuple(sorted((i,(i+1)%B))) for i in range(B)}
 for (a,b),n in ec.items():
  assert n==(1 if (a,b) in boundary else 2)
  if (a,b) not in boundary:assert directed[a,b]==directed[b,a]==1
 for i in range(B):assert directed[i,(i+1)%B]==1 and directed[(i+1)%B,i]==0
 reach={0};stack=[0]
 while stack:
  u=stack.pop()
  for a,b in ec:
   if a==u and b not in reach:reach.add(b);stack.append(b)
   if b==u and a not in reach:reach.add(a);stack.append(a)
 assert reach==set(range(B+N))
 for v,ee in links.items():
  adj=defaultdict(list)
  for a,b in ee:adj[a].append(b);adj[b].append(a)
  seen={next(iter(adj))};stack=list(seen)
  while stack:
   for j in adj[stack.pop()]:
    if j not in seen:seen.add(j);stack.append(j)
  assert len(seen)==len(adj)
  if v<B:
   assert sorted(j for j,ns in adj.items() if len(ns)==1)==sorted([(v-1)%B,(v+1)%B])
   assert all(len(ns) in (1,2) for ns in adj.values())
  else:assert all(len(ns)==2 for ns in adj.values())
 out={'verified':True,'vertices':B+N,'edges':len(ec),'faces':len(triangles),'boundary_caps':B,'exact_equality_center_faces':equality_faces,'strict_center_faces':slack,'minimum_orientation':str(minarea),'minimum_arc_wedge_margin':str(minarc),'oriented_triangles':oriented,'root_sha256':hashlib.sha256((R/'root20.json').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (R/'upper20_verified.json').write_text(json.dumps(out,indent=2));print('UPPER COVER VERIFIED',len(triangles),'triangles,',equality_faces,'equality faces; min orientation',float(minarea),'sec',time.time()-st,flush=True)
 return out
if __name__=='__main__':verify()
