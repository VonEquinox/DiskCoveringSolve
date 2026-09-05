"""Exact upper covering: 36 triangles and ten circular caps.
Six zero-multiplier equality faces are handled by exact symmetry identities,
not by floating-point or straddling interval tests.
"""
from fractions import Fraction as F
from collections import Counter,defaultdict
from pathlib import Path
import json,time
import system14 as S
import symmetry14 as I
from exact_arcs import sincos_turn,SCALE
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def unit(u):
 co,si=sincos_turn(u)
 return [(F(co[0],SCALE),F(co[1],SCALE)),(F(si[0],SCALE),F(si[1],SCALE))]

def combine(*terms):
 d=Counter()
 for fac,v in terms:
  for k,x in v.items():d[k]+=fac*x
 return {k:v for k,v in d.items() if v}
def basis(i):return {i:1}
# All subscripts here are graph nodes: q=0..9, c=10..23, p=24..31.
SUB={9:{10:2,0:-1},15:{10:-1},4:{0:1,10:-2},5:{0:-1},
     25:{0:1,20:1,10:-1},27:{0:1,22:1,10:-1},31:{0:-1,20:1,10:1},29:{0:-1,22:1,10:1}}
def reduce(v):
 out={}
 for i,x in v.items():out=combine((1,out),(x,SUB.get(i,basis(i))))
 return out
def veval(v,p):
 z=[I.point(0),I.point(0)]
 for i,k in v.items():z=I.vadd(z,[I.scale(a,k) for a in p[i]])
 return z

def verify():
 st=time.time();d,z,iv,p=I.rootbox();tlo,thi=iv[S.TID]
 # Symmetry is independently replayed; its proof provides exactly SUB above.
 I.verify()
 edges={frozenset(e) for e in S.EDGES}
 active_w={tuple(f):24+i for i,f in enumerate(S.ACTIVE)}
 formulas={
  (0,1,10):{11:1,10:1,0:-1},(0,9,10):{19:1,10:1,9:-1},
  (4,5,12):{14:1,15:1,4:-1},(5,6,12):{16:1,15:1,5:-1},
  (10,11,12):{21:1,20:1,25:-1},(10,12,13):{23:1,20:1,31:-1}}
 rodforms=[]
 for u,v in S.EDGES:
  f=reduce(combine((1,basis(u)),(-1,basis(v))));rodforms.extend([f,{i:-x for i,x in f.items()}])
 equality_faces=0;slack=[]
 for face in S.FACES:
  if face in active_w:
   assert all(frozenset((10+c,active_w[face])) in edges for c in face);equality_faces+=1
  elif face in formulas:
   w=formulas[face]
   for c in face:assert reduce(combine((1,w),(-1,basis(10+c)))) in rodforms,(face,c)
   equality_faces+=1
  else:
   # A rational circumcenter of the rational midpoint triangle, subsequently
   # rounded to 12 decimal places and verified against the entire root box.
   a,b,c=[[z[2*(10+i-1)+j] for j in (0,1)] for i in face]
   u=[2*(b[j]-a[j]) for j in (0,1)];v=[2*(c[j]-a[j]) for j in (0,1)]
   bu=sum(x*x for x in b)-sum(x*x for x in a);bv=sum(x*x for x in c)-sum(x*x for x in a);de=u[0]*v[1]-u[1]*v[0];assert de
   w=[(bu*v[1]-u[1]*bv)/de,(u[0]*bv-bu*v[0])/de];w=[F(round(q*10**12),10**12) for q in w];wiv=[I.point(q) for q in w]
   mar=min(tlo-I.norm2(I.vsub(p[10+i],wiv))[1] for i in face);assert mar>0
   slack.append({'face':list(face),'witness':list(map(str,w)),'margin':str(mar)})
 assert equality_faces==14 and len(slack)==2
 # Ordered, disjoint rational angular enclosures certify the unit anchors.
 aa=json.load(open(R/'anchor_angles14.json'));cent=list(map(F,aa['centers']));h=F(aa['halfwidth'])
 assert len(cent)==9 and 0<cent[0]-h and cent[-1]+h<1
 assert all(cent[i]+h<cent[i+1]-h for i in range(8))
 for i,m in enumerate(cent,1):
  l,u=m-h,m+h;assert u-l<F(1,4)
  assert I.cross(unit(l),p[i])[0]>0 and I.cross(p[i],unit(u))[0]>0 and I.dot(unit(m),p[i])[0]>0
 # Each boundary disk contains its short unit-circle arc and its convex cap.
 minarc=None
 for i in range(10):
  prev=(i-1)%10;ctr=p[10+i]
  a=I.cross(p[prev],ctr)[0];b=I.cross(ctr,p[i])[0];assert a>0 and b>0 and I.dot(p[prev],p[i])[0]>0
  assert frozenset((prev,10+i)) in edges and frozenset((i,10+i)) in edges
  mar=min(a,b);minarc=mar if minarc is None else min(minarc,mar)
 # The underlying abstract disk, with consistently oriented faces.
 triangles=[]
 for i in range(10):triangles.extend([((i-1)%10,i,10+i),(i,10+(i+1)%10,10+i)])
 triangles.extend(tuple(10+c for c in face) for face in S.FACES)
 directed=Counter();ec=Counter();links=defaultdict(list);minarea=None;oriented=[];collapsed=[]
 for tri in triangles:
  a,b,c=tri;de=I.cross(I.vsub(p[b],p[a]),I.vsub(p[c],p[a]))
  if de[1]<0:b,c=c,b;de=I.cross(I.vsub(p[b],p[a]),I.vsub(p[c],p[a]))
  if tri in [(9,0,10),(4,5,15)]:
   assert not reduce(combine((2,basis(c)),(-1,basis(a)),(-1,basis(b))))
   collapsed.append(list(tri))
  else:
   assert de[0]>0,(tri,de);minarea=de[0] if minarea is None else min(minarea,de[0])
  tri=(a,b,c);oriented.append(list(tri))
  for j in range(3):
   u,v=tri[j],tri[(j+1)%3];directed[u,v]+=1;ec[tuple(sorted((u,v)))]+=1;links[u].append((tri[(j+1)%3],tri[(j+2)%3]))
 assert set(links)==set(range(24)) and len(triangles)==36 and len(ec)==59 and 24-59+36==1
 boundary={tuple(sorted((i,(i+1)%10))) for i in range(10)}
 for (a,b),n in ec.items():
  assert n==(1 if (a,b) in boundary else 2)
  if (a,b) not in boundary:assert directed[a,b]==directed[b,a]==1
 for i in range(10):assert directed[i,(i+1)%10]==1 and directed[(i+1)%10,i]==0
 for v,ee in links.items():
  adj=defaultdict(list)
  for a,b in ee:adj[a].append(b);adj[b].append(a)
  seen={next(iter(adj))};stack=list(seen)
  while stack:
   for j in adj[stack.pop()]:
    if j not in seen:seen.add(j);stack.append(j)
  assert len(seen)==len(adj)
  if v<10:
   assert sorted(j for j,ns in adj.items() if len(ns)==1)==sorted([(v-1)%10,(v+1)%10])
   assert all(len(ns) in (1,2) for ns in adj.values())
  else:assert all(len(ns)==2 for ns in adj.values())
 out={'verified':True,'vertices':24,'edges':59,'faces':36,'boundary_caps':10,'collapsed_ring_triangles':collapsed,'exact_equality_center_faces':14,'strict_center_faces':slack,'minimum_orientation':str(minarea),'minimum_arc_wedge_margin':str(minarc),'oriented_triangles':oriented,'seconds':time.time()-st}
 (R/'upper14_verified.json').write_text(json.dumps(out,indent=2));print('UPPER COVER VERIFIED',len(triangles),'triangles,',equality_faces,'equality faces; sec',time.time()-st,flush=True)
 return out
if __name__=='__main__':verify()
