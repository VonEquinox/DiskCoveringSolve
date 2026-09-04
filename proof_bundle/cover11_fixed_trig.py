#!/usr/bin/env python3
from fractions import Fraction as F
import math
QBITS=160;Q=1<<QBITS

def floordiv(a,b):return a//b
def ceildiv(a,b):return -((-a)//b)
class I:
 __slots__=('lo','hi')
 def __init__(self,lo,hi=None):self.lo=int(lo);self.hi=int(lo if hi is None else hi);assert self.lo<=self.hi
 def __add__(self,o):
  o=o if isinstance(o,I) else I(o);return I(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+(-o)
 def __mul__(self,o):
  o=o if isinstance(o,I) else I(o);p=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
  return I(floordiv(min(p),Q),ceildiv(max(p),Q))
 __rmul__=__mul__
 def divint(self,d):assert d>0;return I(floordiv(self.lo,d),ceildiv(self.hi,d))

def enc(x:F):
 y=x*Q;return I(y.numerator//y.denominator,ceildiv(y.numerator,y.denominator))
def point_nearest(x:F):
 y=x*Q;n=y.numerator//y.denominator;r=y-F(n)
 return n+(1 if r>=F(1,2) else 0)
def upint(x:F):
 y=x*Q;return ceildiv(y.numerator,y.denominator)

# Uniform Lagrange remainders after degree 61/60 on |x|<=7.
SIN_REM=ceildiv(7**62*Q,math.factorial(62))
COS_REM=ceildiv(7**61*Q,math.factorial(61))

def sin_iv(X:I):
 x2=X*X;term=X;s=term
 for k in range(1,31):
  term=(term*x2).divint((2*k)*(2*k+1));s=s-term if k%2 else s+term
 return I(s.lo-SIN_REM,s.hi+SIN_REM)
def cos_iv(X:I):
 x2=X*X;term=I(Q);s=term
 for k in range(1,31):
  term=(term*x2).divint((2*k-1)*(2*k));s=s-term if k%2 else s+term
 return I(s.lo-COS_REM,s.hi+COS_REM)

if __name__=='__main__':
 import mpmath as mp,random
 mp.mp.dps=80
 for z in [0,.1,1,2,3.14,4,6.2]:
  x=F(str(z));iv=enc(x);si=sin_iv(iv);ci=cos_iv(iv)
  sv=mp.sin(mp.mpf(str(z)));cv=mp.cos(mp.mpf(str(z)))
  assert F(si.lo,Q)<=F(str(sv))<=F(si.hi,Q)
  assert F(ci.lo,Q)<=F(str(cv))<=F(ci.hi,Q)
  print(z,float(F(si.hi-si.lo,Q)),float(F(ci.hi-ci.lo,Q)))
