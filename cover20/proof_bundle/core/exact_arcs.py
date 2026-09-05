"""Outward-rounded integer intervals for sin/cos(2*pi*u), u rational.
All decisions use integers and Fraction. No math.sin/cos/pi calls.
"""
from fractions import Fraction as F
if not __debug__:
 raise RuntimeError("Exact checks require Python without -O/PYTHONOPTIMIZE")
from math import factorial
from functools import lru_cache
SCALE=10**60

def ceildiv(a,b):return -((-a)//b)
def floor_fraction(v):return v.numerator//v.denominator
def ceil_fraction(v):return ceildiv(v.numerator,v.denominator)
def arctan_bounds(n,N=100):
 s=sum((F((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(N)),F(0))
 term=F(1,(2*N+1)*n**(2*N+1))
 return (s,s+term) if N%2==0 else (s-term,s)
a,b=arctan_bounds(5);c,d=arctan_bounds(239)
# Machin's identity: pi = 16 atan(1/5) - 4 atan(1/239).
PI=(floor_fraction((16*a-4*d)*SCALE),ceil_fraction((16*b-4*c)*SCALE))
assert 3*SCALE<PI[0]<=PI[1]<4*SCALE

def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
 vals=[x*y for x in a for y in b]
 return min(vals)//SCALE,ceildiv(max(vals),SCALE)
def square(a):
 lo=0 if a[0]<=0<=a[1] else min(a[0]*a[0],a[1]*a[1])
 return lo//SCALE,ceildiv(max(a[0]*a[0],a[1]*a[1]),SCALE)
def divpos(a,n):assert n>0;return a[0]//n,ceildiv(a[1],n)
def rational_interval(v):return floor_fraction(v*SCALE),ceil_fraction(v*SCALE)

@lru_cache(maxsize=50000)
def sincos_turn(u):
 u=F(u)%1;k=floor_fraction(4*u+F(1,2));v=u-F(k,4)
 assert -F(1,8)<=v<=F(1,8)
 xx=[2*F(p,SCALE)*v for p in PI]
 x=(floor_fraction(min(xx)*SCALE),ceil_fraction(max(xx)*SCALE))
 assert max(abs(x[0]),abs(x[1]))<SCALE
 x2=square(x);tc=(SCALE,SCALE);ts=x;co=tc;si=ts
 for j in range(1,20):
  tc=neg(divpos(mul(tc,x2),(2*j-1)*(2*j)))
  ts=neg(divpos(mul(ts,x2),(2*j)*(2*j+1)))
  co=add(co,tc);si=add(si,ts)
 err=ceildiv(SCALE,factorial(40))
 co=(co[0]-err,co[1]+err);si=(si[0]-err,si[1]+err)
 return [(co,si),(neg(si),co),(neg(co),neg(si)),(si,neg(co))][k%4]

QA=10**15;QB=10**18

def nearint(v):
 v=F(v);return (2*v.numerator+v.denominator)//(2*v.denominator)

@lru_cache(maxsize=50000)
def arc_halfplane(lo,hi):
 lo,hi=F(lo),F(hi);assert lo<=hi and hi-lo<F(1,2)
 co,si=sincos_turn((lo+hi)/2)
 a=[nearint(F((v[0]+v[1])*QA,2*SCALE)) for v in (co,si)]
 for ai,iv in zip(a,(co,si)):
  assert (ai-1)*(SCALE//QA)<=iv[0]<=iv[1]<=(ai+1)*(SCALE//QA)
 ch,_=sincos_turn((hi-lo)/2)
 # |a-e(mid)|_infinity <=1/QA, so a.q >=cos(halfwidth)-2/QA.
 b=(ch[0]-2*(SCALE//QA))*QB//SCALE
 return a,b

def dot_interval(u,v):return add(mul(u[0],v[0]),mul(u[1],v[1]))
