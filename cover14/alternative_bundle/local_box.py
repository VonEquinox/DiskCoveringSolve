"""Rigorous anchor-distance bound on an entire angular box (turns)."""
from fractions import Fraction as F
from functools import lru_cache
from exact_arcs import sincos_turn,SCALE
from verify_root import load
from anchor_isolation import RADIUS
@lru_cache(None)
def root_anchor_boxes():
 c,x,Y,rho=load();return [(x[2*i],x[2*i+1],rho) for i in range(9)]
@lru_cache(maxsize=50000)
def trigbox(lo,hi):
 assert lo<=hi and hi-lo<F(1,2)
 pts=[lo,hi]
 for k in range((4*lo).numerator//(4*lo).denominator,(4*hi).numerator//(4*hi).denominator+2):
  p=F(k,4)
  if lo<=p<=hi:pts.append(p)
 vals=[sincos_turn(u) for u in pts]
 return [(min(v[j][0] for v in vals),max(v[j][1] for v in vals)) for j in range(2)]
def bound(lo,hi):
 assert len(lo)==len(hi)==9
 ans=F(0)
 for l,h,(x,y,rho) in zip(lo,hi,root_anchor_boxes()):
  for iv,c in zip(trigbox(F(l),F(h)),(x,y)):
   ans+=max(abs(F(iv[0],SCALE)-(c+rho)),abs(F(iv[1],SCALE)-(c-rho)))**2
 return ans

def accepts(lo,hi):return bound(lo,hi)<=RADIUS**2
