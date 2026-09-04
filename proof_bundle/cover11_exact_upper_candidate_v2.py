#!/usr/bin/env python3
"""Finite rational interval audit that the algebraic full-KKT candidate covers
the unit disk with radius sqrt(t*).

Equal-distance facts are equations of the certified full KKT system.  This
script checks every strict combinatorial/orientation/slack inequality over its
rational root box.
"""
from fractions import Fraction as F
import json
K=json.load(open('/mnt/data/cover11_full_kkt_certificate.json'))
x=[F(int(a),int(b)) for a,b in K['midpoint']];rad=F(int(K['radius'][0]),int(K['radius'][1]));nv={int(k):tuple(v) for k,v in K['nodevar'].items()};edges={tuple(sorted(map(int,e))) for e in K['edges']};tidx=int(K['t_index'])
class IV:
 __slots__=('lo','hi')
 def __init__(self,a,b=None):self.lo=F(a);self.hi=F(a if b is None else b);assert self.lo<=self.hi
 def __add__(self,o):o=o if isinstance(o,IV) else IV(o);return IV(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return IV(-self.hi,-self.lo)
 def __sub__(self,o):return self+(-o)
 def __rsub__(self,o):return IV(o)-self
 def __mul__(self,o):
  o=o if isinstance(o,IV) else IV(o);p=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];return IV(min(p),max(p))
 __rmul__=__mul__
 def sq(self):
  if self.lo<=0<=self.hi:return IV(0,max(self.lo*self.lo,self.hi*self.hi))
  return IV(min(self.lo*self.lo,self.hi*self.hi),max(self.lo*self.lo,self.hi*self.hi))

def point(node):
 if node==4:return (IV(-1),IV(0))
 a,b=nv[node];return (IV(x[a]-rad,x[a]+rad),IV(x[b]-rad,x[b]+rad))
def pmid(node):
 if node==4:return (F(-1),F(0))
 a,b=nv[node];return (x[a],x[b])
P=[point(i) for i in range(29)];M=[pmid(i) for i in range(29)];T=IV(x[tidx]-rad,x[tidx]+rad)
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def orient(a,b,c):return cross(sub(b,a),sub(c,a))
def d2(a,b):q=sub(a,b);return q[0].sq()+q[1].sq()
def msub(a,b):return (a[0]-b[0],a[1]-b[1])
def mcross(a,b):return a[0]*b[1]-a[1]*b[0]
def morient(a,b,c):return mcross(msub(b,a),msub(c,a))

bcent=[10,12,14,16,18,19,17,15,13] # full-node numbers: center c_i is node 9+i
# Check against historical center order [1,3,5,7,9,10,8,6,4].
assert [z-9 for z in bcent]==[1,3,5,7,9,10,8,6,4]
active_cent_faces=[(0,1,3),(0,1,4),(0,2,5),(0,2,6),(0,3,5),(0,4,6),(2,5,7),(2,6,8),(2,9,10)]
slack_cent_faces=[(2,7,9),(2,8,10)]
core=[tuple(9+v for v in f) for f in active_cent_faces+slack_cent_faces]
# Every active face has its three graph edges to the corresponding witness.
for k,f in enumerate(active_cent_faces):
 w=20+k
 for v in f:assert tuple(sorted((w,9+v))) in edges
# Every boundary switch has graph edges to its adjacent boundary centers.
for i in range(9):
 assert tuple(sorted((i,bcent[i]))) in edges
 assert tuple(sorted((i,bcent[(i+1)%9]))) in edges

# Orient all 29 abstract faces counterclockwise by the rational midpoint.
def ccw(f):
 a,b,c=f
 return f if morient(M[a],M[b],M[c])>0 else (a,c,b)
faces=[ccw(f) for f in core]
# Annular triangulation: B_i joins center edge to q_i; A_i joins center i to q_{i-1},q_i.
for i in range(9):
 faces.append(ccw((bcent[i],bcent[(i+1)%9],i)))
 faces.append(ccw((bcent[i],(i-1)%9,i)))
assert len(faces)==29
# Abstract oriented disk check: internal oriented edges cancel; q-cycle is boundary.
inc={}
for a,b,c in faces:
 for e in ((a,b),(b,c),(c,a)):inc[e]=inc.get(e,0)+1
# no same oriented edge repeated
assert max(inc.values())==1
und={}
for a,b,c in faces:
 for u,v in ((a,b),(b,c),(c,a)):
  key=tuple(sorted((u,v)));und.setdefault(key,[]).append((u,v))
for key,arr in und.items():
 if len(arr)==1:
  assert key[0]<9 and key[1]<9 and (key[1]-key[0])%9 in (1,8),key
 elif len(arr)==2:
  assert arr[0]==arr[1][::-1],(key,arr)
 else:raise AssertionError((key,arr))
assert sum(len(v)==1 for v in und.values())==9
# Full combinatorial disk checks: V-E+F=1, connected face adjacency, and
# every vertex link is a cycle (interior) or a path (boundary).
verts={v for f in faces for v in f};assert verts==set(range(20))
assert len(verts)-len(und)+len(faces)==1
face_adj=[set() for _ in faces]
edge_faces={}
for fi,f in enumerate(faces):
 for e in (tuple(sorted((f[0],f[1]))),tuple(sorted((f[1],f[2]))),tuple(sorted((f[2],f[0])))):
  edge_faces.setdefault(e,[]).append(fi)
for arr in edge_faces.values():
 if len(arr)==2:
  a,b=arr;face_adj[a].add(b);face_adj[b].add(a)
seen={0};stack=[0]
while stack:
 u=stack.pop()
 for v in face_adj[u]:
  if v not in seen:seen.add(v);stack.append(v)
assert len(seen)==len(faces)
for v in range(20):
 link={}
 def ladd(a,b):link.setdefault(a,set()).add(b);link.setdefault(b,set()).add(a)
 for f in faces:
  if v in f:
   a,b=[x for x in f if x!=v];ladd(a,b)
 assert link
 lseen={next(iter(link))};lstack=list(lseen)
 while lstack:
  a=lstack.pop()
  for b in link[a]:
   if b not in lseen:lseen.add(b);lstack.append(b)
 assert len(lseen)==len(link)
 deg=sorted(map(len,link.values()));ledges=sum(deg)//2
 if v<9:
  assert deg.count(1)==2 and deg.count(2)==len(deg)-2 and ledges==len(link)-1,(v,deg)
 else:
  assert all(d==2 for d in deg) and ledges==len(link),(v,deg)
# Strict orientation throughout the root box gives an orientation-preserving PL embedding.
omin=None
for f in faces:
 z=orient(P[f[0]],P[f[1]],P[f[2]]);assert z.lo>0,(f,float(z.lo),float(z.hi));omin=z.lo if omin is None or z.lo<omin else omin
# All disk centers lie strictly in the unit disk.
center_margin=None
for node in range(9,20):
 n2=P[node][0].sq()+P[node][1].sq();mar=1-n2.hi;assert mar>0;center_margin=mar if center_margin is None or mar<center_margin else center_margin
# The two nonactive core triangles have circumradius strictly below sqrt(t*).
# Rcirc^2=ab*bc*ca/(4 cross^2).
slack_min=None
for f0 in slack_cent_faces:
 f=tuple(9+v for v in f0);a,b,c=[P[v] for v in f];ab=d2(a,b);bc=d2(b,c);ca=d2(c,a);cr=orient(a,b,c)
 lhs=ab*bc*ca;rhs=4*T*cr.sq();mar=rhs.lo-lhs.hi;assert mar>0,(f,float(mar));slack_min=mar if slack_min is None or mar<slack_min else slack_min
# The nine boundary anchors themselves occur in the claimed cyclic order.
# Since every consecutive cross product is positive and each gap is already
# constrained to be below pi by the adjacent-dot tests below, the polygonal
# boundary q0,...,q8 is the positively oriented simple convex 9-gon inscribed
# in the unit circle.
qcyclemin=None
for i in range(9):
 z=cross(P[i],P[(i+1)%9]);assert z.lo>0,(i,float(z.lo),float(z.hi))
 qcyclemin=z.lo if qcyclemin is None or z.lo<qcyclemin else qcyclemin
# Explicit half-plane signs remove any possible ambiguity in the wrap at q4.
for i in (0,1,7,8):assert P[i][0].lo>0
for i in (2,3,5,6):assert P[i][0].hi<0
for i in (0,1,2,3):assert P[i][1].lo>0
for i in (5,6,7,8):assert P[i][1].hi<0
assert P[4][0].lo==P[4][0].hi==-1 and P[4][1].lo==P[4][1].hi==0
# Correct cyclic minor arcs and center direction: endpoints q_{i-1},q_i,
# center b_i.  These strict signs imply the center direction lies inside the
# minor arc, and dot>0 makes the arc shorter than pi/2 (hence <pi).
arcmin=None
for i,cnode in enumerate(bcent):
 qp=P[(i-1)%9];qn=P[i];c=P[cnode]
 vals=[cross(qp,c),cross(c,qn),dot(qp,qn)]
 for z in vals:
  assert z.lo>0,(i,float(z.lo),float(z.hi));arcmin=z.lo if arcmin is None or z.lo<arcmin else arcmin
# Exact algebraic candidate lies below the previous rational upper radius.
Rrat=F(379983853121,10**12);assert T.hi<Rrat*Rrat
print('EXACT ALGEBRAIC CANDIDATE COVER VERIFIED')
print('disk complex V/E/F',len(verts),len(und),len(faces))
print('faces',len(faces),'min orientation',float(omin),'center-unit margin',float(center_margin))
print('slack circum polynomial margin',float(slack_min),'q-cycle margin',float(qcyclemin),'arc sign margin',float(arcmin))
print('t interval',float(T.lo),float(T.hi),'old rational R2 margin',float(Rrat*Rrat-T.hi))
