"""Exact 19-disc cover: t=1/13, 24 equilateral center faces, 48 total triangles.
Coordinates are rational in the orthogonal basis (1,0),(0,sqrt(3)).
Degenerate ring triangles are accepted only as genuine exact equalities.
"""
from fractions import Fraction as F
from collections import Counter,defaultdict
from pathlib import Path
import json,time,hashlib
import geometry19 as S
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def quadrant(p):
 x,y=p;assert x or y
 if x>0 and y>=0:return 0
 if y>0 and x<=0:return 1
 if x<0 and y<=0:return 2
 return 3

def verify():
 st=time.time();S.check_geometry();B=S.B;N=S.N;p=S.ANCHORS+S.CENTERS
 qclass=list(map(quadrant,S.ANCHORS))
 assert qclass==sorted(qclass) and qclass[0]==0 and qclass[-1]==3
 assert S.ANCHORS[0]==(1,0)
 for i in range(B):
  a,b=S.ANCHORS[i],S.ANCHORS[(i+1)%B]
  assert S.cross(a,b)>0 and S.dot(a,b)>0
 for i in range(B):
  a,b=S.ANCHORS[(i-1)%B],S.ANCHORS[i];c=S.CENTERS[i]
  assert S.cross(a,c)>0 and S.cross(c,b)>0
  assert S.norm2(S.sub(a,c))==S.norm2(S.sub(b,c))==S.T
 triangles=[]
 for i in range(B):triangles.extend([((i-1)%B,i,B+i),(i,B+(i+1)%B,B+i)])
 for face in S.FACES:triangles.append(tuple(B+i for i in face))
 directed=Counter();ec=Counter();links=defaultdict(list);zero=[];minpos=None;oriented=[]
 for k,tri in enumerate(triangles):
  a,b,c=tri;de=S.cross(S.sub(p[b],p[a]),S.sub(p[c],p[a]))
  if de<0:
   assert k>=2*B # ring orientation is fixed combinatorially.
   b,c=c,b;de=-de
  if de==0:
   assert k<2*B and k%2==0
   i=k//2;assert S.CENTERS[i]==S.scale(S.add(S.ANCHORS[(i-1)%B],S.ANCHORS[i]),F(1,2))
   zero.append(k)
  else:minpos=de if minpos is None else min(minpos,de)
  tri=(a,b,c);oriented.append(list(tri))
  for j in range(3):
   u,v=tri[j],tri[(j+1)%3];directed[u,v]+=1;ec[tuple(sorted((u,v)))]+=1;links[u].append((tri[(j+1)%3],tri[(j+2)%3]))
 assert zero==[0,4,8,12,16,20]
 assert len(triangles)==48 and len(ec)==78 and B+N==31 and B+N-len(ec)+len(triangles)==1
 boundary={tuple(sorted((i,(i+1)%B))) for i in range(B)}
 for (a,b),cnt in ec.items():
  assert cnt==(1 if (a,b) in boundary else 2)
  if (a,b) not in boundary:assert directed[a,b]==directed[b,a]==1
 for i in range(B):assert directed[i,(i+1)%B]==1 and directed[(i+1)%B,i]==0
 reach={0};todo=[0]
 while todo:
  u=todo.pop()
  for a,b in ec:
   if a==u and b not in reach:reach.add(b);todo.append(b)
   if b==u and a not in reach:reach.add(a);todo.append(a)
 assert reach==set(range(B+N)) and set(links)==reach
 for v,ee in links.items():
  adj=defaultdict(list)
  for a,b in ee:adj[a].append(b);adj[b].append(a)
  seen={next(iter(adj))};todo=list(seen)
  while todo:
   for j in adj[todo.pop()]:
    if j not in seen:seen.add(j);todo.append(j)
  assert len(seen)==len(adj)
  if v<B:
   assert sorted(j for j,ns in adj.items() if len(ns)==1)==sorted([(v-1)%B,(v+1)%B])
   assert all(len(ns) in (1,2) for ns in adj.values())
  else:assert all(len(ns)==2 for ns in adj.values())
 # Exact sixty-place radical intervals, certified by integer-square comparisons.
 den=10**60
 from math import isqrt
 a=isqrt(den*den//13);b=isqrt(13*den*den)
 assert 13*a*a<den*den<13*(a+1)*(a+1)
 assert b*b<13*den*den<(b+1)*(b+1)
 def dec(v):return str(v//den)+'.'+str(v%den).zfill(60)
 out={'verified':True,'exact_radius':'1/sqrt(13)','exact_large_radius':'sqrt(13)','t':'1/13','decimal_intervals':{'r_lower':dec(a),'r_upper':dec(a+1),'R_lower':dec(b),'R_upper':dec(b+1)},'vertices':B+N,'edges':len(ec),'triangles':len(triangles),'caps':B,'center_equilateral_faces':24,'degenerate_ring_faces':zero,'minimum_positive_orientation_over_sqrt3':str(minpos),'oriented_triangles':oriented,'geometry_sha256':hashlib.sha256((R/'geometry19.py').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (R/'upper19_verified.json').write_text(json.dumps(out,indent=2));print('EXACT UPPER COVER VERIFIED',out['decimal_intervals'],flush=True);return out
if __name__=='__main__':verify()
