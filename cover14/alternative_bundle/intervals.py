"""Elementary rational interval arithmetic."""
from dataclasses import dataclass
from fractions import Fraction as F
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
