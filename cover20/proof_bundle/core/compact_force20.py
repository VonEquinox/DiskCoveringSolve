"""Integer force certificates in cycle coordinates. Standard library only."""
from fractions import Fraction as F
from math import isqrt
from exact_arcs import QA,QB,arc_halfplane
from functools import lru_cache
if not __debug__:raise RuntimeError('Run without -O')
N=20;V=58;RU=F('0.27084817843999')
def cr(n):
 q=isqrt(n);return q+(q*q<n)
@lru_cache(maxsize=128)
def tree(ed):
 adj=[[] for _ in range(V)]
 for k,(u,v) in enumerate(ed):
  assert 0<=u<V and 0<=v<V and u!=v
  adj[u].append((v,k));adj[v].append((u,k))
 parent=[None]*V;parent[0]=(-1,-1);order=[0];used=set()
 for u in order:
  for v,k in adj[u]:
   if parent[v] is None:parent[v]=(u,k);used.add(k);order.append(v)
 assert len(order)==V
 return parent,order,tuple(k for k in range(len(ed)) if k not in used)
def reconstruct(B,ed,lo,hi,c):
 assert set(c)=={'q','l','b','f'}
 Q=c['q'];la=c['l'];beta=c['b'];cy=c['f']
 assert type(Q)is int and Q>0 and Q<=10**18
 assert len(la)==B-1 and len(beta)==2*(B-1)
 assert all(type(x)is int and x>=0 for x in la)
 assert all(type(x)is int for x in beta+cy)
 parent,order,chord=tree(tuple(ed));assert len(cy)==2*len(chord)
 aa=[];bb=[]
 for l,h in zip(lo,hi):a,b=arc_halfplane(l,h);aa.append(a);bb.append(b)
 assert len(aa)==B-1 and len(lo)==len(hi)==B-1
 balance=[[0,0] for _ in range(V)]
 for i in range(B-1):
  for d in range(2):
   q=la[i]*aa[i][d]-beta[2*i+d]*QA
   balance[i+1][d]=q;balance[0][d]-=q
 flow=[[0,0] for _ in ed]
 for j,k in enumerate(chord):
  u,v=ed[k]
  for d in range(2):
   q=cy[2*j+d]*QA;flow[k][d]=q;balance[u][d]-=q;balance[v][d]+=q
 for v in reversed(order[1:]):
  u,k=parent[v];sg=1 if ed[k][0]==v else -1
  for d in range(2):flow[k][d]=sg*balance[v][d];balance[u][d]+=balance[v][d];balance[v][d]=0
 assert all(z==[0,0] for z in balance)
 # A fresh divergence recomputation checks the reconstruction rather than trusting it.
 div=[[0,0] for _ in range(V)]
 for (u,v),z in zip(ed,flow):
  for d in range(2):div[u][d]+=z[d];div[v][d]-=z[d]
 for i in range(1,V):
  if i<B:assert div[i]==[la[i-1]*aa[i-1][d]-beta[2*(i-1)+d]*QA for d in range(2)]
  else:assert div[i]==[0,0]
 assert div[0]==[-sum(div[i][d] for i in range(1,V)) for d in range(2)]
 return aa,bb,flow

def check(B,ed,lo,hi,c):
 aa,bb,flow=reconstruct(B,ed,lo,hi,c);la=c['l'];beta=c['b'];ns=sum(cr(x*x+y*y) for x,y in flow)
 assert ns>0
 bn=0
 for i in range(B-1):
  x,y=beta[2*i:2*i+2];bn+=la[i]*(bb[i]-aa[i][0]*(QB//QA))-(cr((x*QA)**2+(y*QA)**2)-x*QA)*(QB//QA)
 return F(bn*QA,QB*ns)-RU

def rationalize(B,ed,lo,hi,lam,support,flow):
 parent,order,chords=tree(tuple(ed))
 last=None
 for Q in [10**6,10**8,10**10,10**12,10**15]:
  c={'q':Q,'l':[max(0,int(round(float(x)*Q))) for x in lam],
     'b':[int(round(float(x)*Q)) for row in support for x in row],
     'f':[int(round(float(flow[k][d])*Q)) for k in chords for d in range(2)]}
  try:mar=check(B,ed,lo,hi,c)
  except (AssertionError,ValueError):continue
  last=mar
  if mar>0:return c,mar
 return None,last
