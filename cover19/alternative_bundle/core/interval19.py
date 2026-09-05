"""Rational Cartesian interval arithmetic for exact Q(sqrt(3)) anchors."""
from fractions import Fraction as F
import geometry19 as S
if not __debug__:raise RuntimeError('Run without -O')
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def point(x):return F(x),F(x)
def dot(a,b):return add(mul(a[0],b[0]),mul(a[1],b[1]))
def cross(a,b):return sub(mul(a[0],b[1]),mul(a[1],b[0]))
def vsub(a,b):return [sub(x,y) for x,y in zip(a,b)]
def norm2(a):return dot(a,a)
