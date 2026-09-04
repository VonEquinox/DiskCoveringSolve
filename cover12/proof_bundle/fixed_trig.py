#!/usr/bin/env python3
"""Fast rigorous fixed-point Taylor bounds for sin/cos and affine minorants."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
BITS=190;Q=1<<BITS

def floorF(x):x=F(x)*Q;return x.numerator//x.denominator
def ceilF(x):x=F(x)*Q;return -((-x.numerator)//x.denominator)
def cdiv(a,n):return -((-a)//n)
class Iv:
 __slots__=('l','u')
 def __init__(self,l,u=None):self.l=int(l);self.u=int(l if u is None else u);assert self.l<=self.u
 @staticmethod
 def frac(l,u=None):return Iv(floorF(l),ceilF(l if u is None else u))
 def __add__(self,o):o=iv(o);return Iv(self.l+o.l,self.u+o.u)
 __radd__=__add__
 def __neg__(self):return Iv(-self.u,-self.l)
 def __sub__(self,o):return self+(-iv(o))
 def __mul__(self,o):
  o=iv(o);p=[self.l*o.l,self.l*o.u,self.u*o.l,self.u*o.u];return Iv(min(p)//Q,cdiv(max(p),Q))
 __rmul__=__mul__
 def divint(self,n):assert n>0;return Iv(self.l//n,cdiv(self.u,n))
 def lo(self):return F(self.l,Q)
 def hi(self):return F(self.u,Q)
def iv(x):return x if isinstance(x,Iv) else Iv.frac(x)

def powiv(x,n):
 z=Iv(Q)
 for _ in range(n):z=z*x
 return z
@lru_cache(maxsize=500000)
def sincos_point(x,n=42):
 x=F(x);X=Iv.frac(x);x2=X*X
 ts=X;ss=Iv(0)
 for k in range(n+1):
  ss=ss+ts
  ts=-(ts*x2).divint((2*k+2)*(2*k+3))
 tc=Iv(Q);cc=Iv(0)
 for k in range(n+1):
  cc=cc+tc
  tc=-(tc*x2).divint((2*k+1)*(2*k+2))
 # Taylor remainder theorem (add fixed outward bound)
 ax=abs(x);rs=F(ax**(2*n+3),factorial(2*n+3));rc=F(ax**(2*n+2),factorial(2*n+2))
 rsn=ceilF(rs);rcn=ceilF(rc)
 return F(ss.l-rsn,Q),F(ss.u+rsn,Q),F(cc.l-rcn,Q),F(cc.u+rcn,Q)

def cos_lower_interval(l,u,PI_L,PI_U):
 l,u=F(l),F(u)
 if l<=PI_U and u>=PI_L:return F(-1)
 if u<PI_L:return sincos_point(u)[2]
 if l>PI_U:return sincos_point(l)[2]
 return F(-1)
@lru_cache(maxsize=500000)
def line_minorant(l,u,PI_L,PI_U):
 l,u=F(l),F(u);assert l<=u
 a=(l+u)/2;h=(u-l)/2;sl,su,cl,cu=sincos_point(a);A=(sl+su)/2;eps=(su-sl)/2;phil=1-cu;m=cos_lower_interval(l,u,PI_L,PI_U);B=phil-A*a-eps*h+min(F(0),m)*h*h/2
 return A,B
