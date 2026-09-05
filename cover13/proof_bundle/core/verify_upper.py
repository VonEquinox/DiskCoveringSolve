#!/usr/bin/env python3
"""Exact interval verification that the isolated n=13 KKT root covers the unit disk."""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
from collections import Counter,defaultdict,deque
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import cover13_kkt_system as K
@dataclass(frozen=True)
class I:
 lo:F;hi:F
 def __post_init__(self):assert self.lo<=self.hi
 @staticmethod
 def point(x):x=x if isinstance(x,F) else F(x);return I(x,x)
 def __add__(self,o):o=asI(o);return I(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+(-asI(o))
 def __rsub__(self,o):return asI(o)-self
 def __mul__(self,o):
  o=asI(o);z=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi);return I(min(z),max(z))
 __rmul__=__mul__
 def sq(self):
  if self.lo<=0<=self.hi:return I(F(0),max(self.lo*self.lo,self.hi*self.hi))
  z=(self.lo*self.lo,self.hi*self.hi);return I(min(z),max(z))
 def mid(self):return (self.lo+self.hi)/2
 def width(self):return self.hi-self.lo
def asI(x):return x if isinstance(x,I) else I.point(x)
def pt(x,y):return (asI(x),asI(y))
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def orient(a,b,c):return cross(sub(b,a),sub(c,a))
def dist2(a,b):d=sub(a,b);return d[0].sq()+d[1].sq()
CERT=json.load(open(ROOT/'cover13_krawczyk_cert.json'));ROOTDATA=json.load(open(ROOT/'cover13_kkt_root_119d.json'))
Qx=int(CERT['Qx']);rho=F(int(CERT['rho_num']),int(CERT['rho_den']))
xmid=[F(int(n),Qx) for n in CERT['xnum']];X=[I(x-rho,x+rho) for x in xmid];assert len(X)==K.N==119
def coord(node,a):
 if node==K.FIXED:return I.point(1 if a==0 else 0)
 if node in K.QIDS:return X[K.QSTART+2*K.QIDS.index(node)+a]
 assert node in K.FREEIDS
 return X[K.FSTART+2*K.FREEIDS.index(node)+a]
def pos(node):return (coord(node,0),coord(node,1))
Q=[pos(i) for i in range(10)];C=[pos(10+i) for i in range(13)];P=[pos(23+i) for i in range(9)];t=X[K.TID]
EDGES={tuple(sorted(e)) for e in K.EDGES}
for i in range(10):
 assert tuple(sorted((i,10+i))) in EDGES
 assert tuple(sorted((i,10+(i+1)%10))) in EDGES
ACTIVE=[(0,10,12),(1,2,10),(3,4,11),(3,10,11),(4,5,11),(5,6,11),(6,7,11),(7,11,12),(8,9,12)]
assert [tuple(x) for x in ROOTDATA['active_faces']]==ACTIVE
for k,fa in enumerate(ACTIVE):
 for c in fa:assert tuple(sorted((23+k,10+c))) in EDGES
assert K.QIDS==list(range(1,10))
# Boundary order: each short step is positively oriented and <pi/2. Sign pattern proves one winding.
quadrants=[(1,0),(1,1),(1,1),(-1,1),(-1,1),(-1,-1),(-1,-1),(-1,-1),(1,-1),(1,-1)]
boundary_order_margins=[];center_cone_margins=[]
for i,(sx,sy) in enumerate(quadrants):
 if sx>0:assert Q[i][0].lo>0
 else:assert Q[i][0].hi<0
 if sy>0:assert Q[i][1].lo>0
 elif sy<0:assert Q[i][1].hi<0
 else:assert Q[i][1].lo==Q[i][1].hi==0
for i in range(10):
 j=(i+1)%10;cp=cross(Q[i],Q[j]);dp=dot(Q[i],Q[j]);assert cp.lo>0 and dp.lo>0
 boundary_order_margins += [cp.lo,dp.lo]
for i in range(10):
 p=(i-1)%10;z1=cross(Q[p],C[i]);z2=cross(C[i],Q[i]);assert z1.lo>0 and z2.lo>0
 center_cone_margins += [z1.lo,z2.lo]
ALL=[(0,1,10),(0,9,12),(0,10,12),(1,2,10),(2,3,10),(3,4,11),(3,10,11),(4,5,11),(5,6,11),(6,7,11),(7,8,12),(7,11,12),(8,9,12),(10,11,12)]
assert [tuple(x) for x in ROOTDATA['all_faces']]==ALL
faces=[]
for i in range(10):
 faces.append(((i-1)%10,i,10+i));faces.append((i,10+(i+1)%10,10+i))
for fa in ALL:
 ids=[10+x for x in fa];val=orient(C[fa[0]],C[fa[1]],C[fa[2]])
 if val.hi<0:ids[1],ids[2]=ids[2],ids[1]
 else:assert val.lo>0
 faces.append(tuple(ids))
assert len(faces)==34
def V(v):return Q[v] if v<10 else C[v-10]
orientation_margins=[]
for fa in faces:
 o=orient(V(fa[0]),V(fa[1]),V(fa[2]));assert o.lo>0;orientation_margins.append(o.lo)
ec=Counter()
for fa in faces:
 for a,b in ((fa[0],fa[1]),(fa[1],fa[2]),(fa[2],fa[0])):ec[tuple(sorted((a,b)))]+=1
outer={tuple(sorted((i,(i+1)%10))) for i in range(10)}
assert {e for e,n in ec.items() if n==1}==outer
assert all(n in (1,2) for n in ec.values())
assert len(ec)==56 and 23-len(ec)+len(faces)==1
adj=defaultdict(set)
for a,b in ec:adj[a].add(b);adj[b].add(a)
seen={0};dq=deque([0])
while dq:
 a=dq.popleft()
 for b in adj[a]:
  if b not in seen:seen.add(b);dq.append(b)
assert seen==set(range(23))
# Combinatorial closed-disk audit.  Together with positive face orientations and
# the boundary homeomorphism, the PL degree lemma proves that the affine map is
# a homeomorphism onto the convex anchor decagon; no pairwise noncrossing test is
# needed (and the exact root has harmless collinear extensions on its symmetry axis).
link_report={}
for v in range(23):
 ledges=[]
 for fa in faces:
  if v in fa:
   u=[x for x in fa if x!=v];assert len(u)==2;ledges.append(tuple(sorted(u)))
 ladj=defaultdict(set)
 for a,b in ledges:ladj[a].add(b);ladj[b].add(a)
 assert ladj
 lseen={next(iter(ladj))};ldq=deque(lseen)
 while ldq:
  a=ldq.popleft()
  for b in ladj[a]:
   if b not in lseen:lseen.add(b);ldq.append(b)
 assert lseen==set(ladj)
 deg=sorted(len(z) for z in ladj.values())
 if v<10:
  assert deg.count(1)==2 and all(x in (1,2) for x in deg)
  link_report[str(v)]='path'
 else:
  assert all(x==2 for x in deg)
  link_report[str(v)]='cycle'
INACTIVE_WITNESS={
 (0,1,10):('0.56878191710999470576','0.06679623442179143111'),
 (0,9,12):('0.43033539731774933523','-0.37787041691514644537'),
 (2,3,10):('0.22304103736749614217','0.50521772187193447756'),
 (7,8,12):('-0.10317648754583604032','-0.54253743028581913332'),
 (10,11,12):('0.03404965701116063198','-0.01060132685718953756'),
}
inactive_margins=[];active_index={fa:k for k,fa in enumerate(ACTIVE)}
assert set(ALL)==set(ACTIVE)|set(INACTIVE_WITNESS)
for fa in ALL:
 if fa in active_index:
  k=active_index[fa]
  for c in fa:assert tuple(sorted((23+k,10+c))) in EDGES
 else:
  wx,wy=INACTIVE_WITNESS[fa];w=pt(F(wx),F(wy))
  for c in fa:
   mar=t.lo-dist2(w,C[c]).hi;assert mar>0;inactive_margins.append(mar)
all_strict=boundary_order_margins+center_cone_margins+orientation_margins+inactive_margins
assert all(z>0 for z in all_strict)
out={'verified':True,'root_box_radius':str(rho),'t_interval':[str(t.lo),str(t.hi)],
 'vertices':23,'straight_faces':34,'straight_edges':56,'boundary_faces':20,'center_faces':14,
 'active_center_faces':9,'inactive_center_faces':5,
 'min_boundary_order_margin':float(min(boundary_order_margins)),
 'min_center_cone_margin':float(min(center_cone_margins)),
 'min_triangle_orientation_margin':float(min(orientation_margins)),
 'min_inactive_witness_margin':float(min(inactive_margins)),
 'krawczyk_cert_sha256':hashlib.sha256((ROOT/'cover13_krawczyk_cert.json').read_bytes()).hexdigest(),
 'root_data_sha256':hashlib.sha256((ROOT/'cover13_kkt_root_119d.json').read_bytes()).hexdigest()}
(ROOT/'cover13_upper_verified.json').write_text(json.dumps(out,indent=2))
print('COVER13 UPPER COVER VERIFIED');print(json.dumps(out,indent=2))
