"""Exact anchor-only local leaf, for the 42-rod rational candidate."""
from fractions import Fraction as F
from functools import lru_cache
import geometry19 as S
import interval19 as I
from exact_arcs import sincos_turn,SCALE
from anchor_isolation19 import RADIUS
if not __debug__:raise RuntimeError('Run without -O')
@lru_cache(maxsize=1)
def _anchors():return S.cartesian_intervals(S.ANCHORS)
def local(lo,hi):
 assert len(lo)==len(hi)==S.B-1
 anchors=_anchors();ds=F(0)
 for i,(l,h) in enumerate(zip(lo,hi),1):
  if h-l>=F(1,2):return None
  def q(u):
   co,si=sincos_turn(u);return [(F(a,SCALE),F(b,SCALE)) for a,b in (co,si)]
  if I.dot(q((l+h)/2),anchors[i])[0]<=0:return None
  low=min(I.dot(q(l),anchors[i])[0],I.dot(q(h),anchors[i])[0]);ds+=2-2*low
 return RADIUS*RADIUS-ds
