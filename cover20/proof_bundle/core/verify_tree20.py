"""Exact complete-tree acceptance, including every incoming angular box."""
from fractions import Fraction as F
import domain20 as D
import compact_force20 as C
from candidate_global20 import local
if not __debug__:raise RuntimeError('Run without -O')
def verify(B,faces,nodes,candidate=False):
 assert type(nodes)is list and nodes
 mat=D.bounds(B,faces);lo,hi=D.rootbox(mat);ed=D.edges(B,faces)
 seen=bytearray(len(nodes));stack=[(0,lo,hi,0)];counts={'SPLIT':0,'DUAL':0,'LOCAL':0,'EMPTY':0};mind=minl=None;depthmax=0;scales={}
 while stack:
  k,lo,hi,dep=stack.pop();assert type(k)is int and 0<=k<len(nodes) and not seen[k],'invalid, repeated or cyclic node'
  seen[k]=1;n=nodes[k];assert type(n)is dict and n.get('kind') in counts
  kind=n['kind'];counts[kind]+=1;depthmax=max(depthmax,dep);box=D.tighten(mat,lo,hi)
  if kind=='EMPTY':assert set(n)=={'kind'} and box is None;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   assert set(n)=={'kind','axis','mid','children'};axis=n['axis'];ch=n['children'];assert type(axis)is int and 0<=axis<B-1 and type(ch)is list and len(ch)==2 and ch[0]!=ch[1]
   assert all(type(v)is int and k<v<len(nodes) for v in ch)
   mid=F(n['mid']);assert lo[axis]<mid<hi[axis]
   h=list(hi);h[axis]=mid;l=list(lo);l[axis]=mid
   stack.extend([(ch[1],l,hi,dep+1),(ch[0],lo,h,dep+1)])
  elif kind=='DUAL':
   assert set(n)=={'kind','cert'};mar=C.check(B,ed,lo,hi,n['cert']);assert mar>0,('nonpositive force margin',k,mar)
   mind=mar if mind is None else min(mind,mar);q=str(n['cert']['q']);scales[q]=scales.get(q,0)+1
  elif kind=='LOCAL':
   assert set(n)=={'kind'} and candidate and B==13
   mar=local(tuple(lo),tuple(hi));assert mar is not None and mar>0
   minl=mar if minl is None else min(minl,mar)
 assert all(seen),'unreachable nodes'
 assert counts['SPLIT']+1==counts['DUAL']+counts['LOCAL']+counts['EMPTY']
 return {'nodes':len(nodes),'counts':counts,'maxdepth':depthmax,'minimum_force_margin':mind,'minimum_local_margin':minl,'scales':scales}
