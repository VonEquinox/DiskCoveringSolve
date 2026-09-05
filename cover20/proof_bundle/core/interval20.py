"""Rational interval operations and the certified n=20 root box."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
import system20 as S
if not __debug__:raise RuntimeError('Run without -O')
R=S.R

def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[a0*b0 for a0 in a for b0 in b];return min(v),max(v)
def point(x):x=F(x);return x,x
def scale(a,x):return mul(a,point(x))
def vadd(a,b):return [add(x,y) for x,y in zip(a,b)]
def vsub(a,b):return [sub(x,y) for x,y in zip(a,b)]
def dot(a,b):return add(mul(a[0],b[0]),mul(a[1],b[1]))
def cross(a,b):return sub(mul(a[0],b[1]),mul(a[1],b[0]))
def norm2(a):return dot(a,a)

def rootbox(rho=F(1,10**80)):
 d=json.load(open(R/'root20.json'));q=int(d['Qx']);z=[F(int(v),q) for v in d['xnum']];iv=[(v-rho,v+rho) for v in z];p=[[point(1),point(0)]]+[iv[2*i:2*i+2] for i in range(S.G//2)]
 return d,z,iv,p
