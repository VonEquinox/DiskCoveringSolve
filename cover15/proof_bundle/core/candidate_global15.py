"""Candidate topology: exact force exclusions plus a 13/100 anchor-only theorem."""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import json,gzip,time,hashlib
import system15 as S
import interval15 as I
from exact_arcs import sincos_turn,SCALE
R=S.R
if not __debug__:raise RuntimeError('Run without -O')

@lru_cache(maxsize=1)
def _anchors():return I.rootbox()[3]

def local(lo,hi):
 # A point on the unit circle has squared chord distance 2-2*dot.
 # On a subsemicircle whose midpoint dot is positive, the minimum is at an
 # endpoint (the only interior dot minimum is -1, which is excluded).
 p=_anchors();ds=F(0)
 for i,(l,h) in enumerate(zip(lo,hi),1):
  if h-l>=F(1,2):return None
  def q(u):
   co,si=sincos_turn(u);return [(F(a,SCALE),F(b,SCALE)) for a,b in (co,si)]
  if I.dot(q((l+h)/2),p[i])[0]<=0:return None
  dotmin=min(I.dot(q(l),p[i])[0],I.dot(q(h),p[i])[0]);ds+=2-2*dotmin
 return F(169,10000)-ds
