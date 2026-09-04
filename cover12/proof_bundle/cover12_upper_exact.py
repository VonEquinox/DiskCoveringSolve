#!/usr/bin/env python3
"""Exact interval verification of the D3 candidate upper covering geometry."""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from exact_core import PI_L,PI_U
from fixed_trig import sincos_point
BITS=230;Q=1<<BITS

def ceildiv(a,b):return -((-a)//b)
def ffloor(x):z=F(x)*Q;return z.numerator//z.denominator
def fceil(x):z=F(x)*Q;return ceildiv(z.numerator,z.denominator)
class Iv:
 __slots__=('l','u')
 def __init__(self,l,u=None):self.l=int(l);self.u=int(l if u is None else u);assert self.l<=self.u
 @staticmethod
 def frac(l,u=None):return Iv(ffloor(l),fceil(l if u is None else u))
 def __add__(self,o):o=iv(o);return Iv(self.l+o.l,self.u+o.u)
 __radd__=__add__
 def __neg__(self):return Iv(-self.u,-self.l)
 def __sub__(self,o):return self+(-iv(o))
 def __rsub__(self,o):return iv(o)-self
 def __mul__(self,o):
  o=iv(o);p=[self.l*o.l,self.l*o.u,self.u*o.l,self.u*o.u];return Iv(min(p)//Q,ceildiv(max(p),Q))
 __rmul__=__mul__
 def inv(self):
  assert self.l>0;return Iv((Q*Q)//self.u,ceildiv(Q*Q,self.l))
 def __truediv__(self,o):return self*iv(o).inv()
 def lo(self):return F(self.l,Q)
 def hi(self):return F(self.u,Q)
 def width(self):return F(self.u-self.l,Q)
def iv(x):return x if isinstance(x,Iv) else Iv.frac(x)
ZERO=Iv(0);HALF=Iv.frac(F(1,2))

def load():
 d=json.load(open(ROOT/'cover12_symmetric_kkt_cert.json'));q=int(d['Qx']);rho=F(int(d['rho_num']),int(d['rho_den']));m={n:F(int(a),q) for n,a in zip(d['names'],d['xnum'])};X={n:Iv.frac(a-rho,a+rho) for n,a in m.items()};return X,m,rho

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def mul(s,a):return (s*a[0],s*a[1])
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def orient(a,b,c):return cross(sub(b,a),sub(c,a))
def rot120(a,z):return (-(a[0]+z*a[1])*HALF,(z*a[0]-a[1])*HALF)
def rot240(a,z):return ((-a[0]+z*a[1])*HALF,-(z*a[0]+a[1])*HALF)
def reflect60(a,z):return ((-a[0]+z*a[1])*HALF,(z*a[0]+a[1])*HALF)

def circum2(a,b,c):
 ba=sub(b,a);ca=sub(c,a);db=dot(b,b)-dot(a,a);dc=dot(c,c)-dot(a,a);den=2*cross(ba,ca);assert den.l>0
 ux=(db*ca[1]-dc*ba[1])/den;uy=(ba[0]*dc-ca[0]*db)/den;return dot(sub((ux,uy),a),sub((ux,uy),a))

def geometry(X):
 r,u,v,w,c,s,z=[X[k] for k in ['r','u','v','w','c','s','z']]
 gp=(u*c,u*s);gm=(u*c,-u*s)
 rawC=[gp,gm,rot120(gp,z),rot120(gm,z),rot240(gp,z),rot240(gm,z),(-v,ZERO),rot120((-v,ZERO),z),rot240((-v,ZERO),z),(w,ZERO),rot120((w,ZERO),z),rot240((w,ZERO),z)]
 C=[rawC[i] for i in [6,5,4,7,1,0,8,3,2,9,10,11]]
 c2=c*c-s*s;s2=2*c*s;q0=(iv(1),ZERO);q1=(c2,s2);q2=reflect60(q1,z)
 rawQ=[q0,q1,q2,rot120(q0,z),rot120(q1,z),rot120(q2,z),rot240(q0,z),rot240(q1,z),rot240(q2,z)]
 Qp=[rawQ[i] for i in [5,6,7,8,0,1,2,3,4]]
 return C,Qp

def verify():
 t=time.time();X,m,rho=load();C,Qp=geometry(X);r=X['r'];c=X['c'];s=X['s'];z=X['z']
 assert r.hi()<F(1,2) and max(X[k].hi() for k in ['u','v','w'])<1 and X['w'].hi()<r.lo()
 # h and k=pi/3-2h lie in (0,pi/2), and both corresponding half-gaps < DELTA/2.
 Delta=F(73891,100000);sinDl=sincos_point(Delta/2)[0];L=(c*c-s*s)*HALF+z*c*s;M=z*(c*c-s*s)*HALF-c*s
 assert min(c.lo(),s.lo(),L.lo(),M.lo())>0
 assert s.hi()<sinDl and M.hi()<sinDl
 # 31-face oriented PL disk
 core=[(0,1,11),(0,10,8),(0,11,10),(1,2,11),(2,3,11),(3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),(6,10,9),(7,8,10),(9,10,11)]
 allfaces=list(core)+[(12+(i-1)%9,12+i,i) for i in range(9)]+[(12+i,(i+1)%9,i) for i in range(9)]
 assert len(allfaces)==31
 # Purely combinatorial verification of the oriented disk triangulation.
 directed={};undirected={}
 for f in allfaces:
  assert len(set(f))==3
  for aa,bb in ((f[0],f[1]),(f[1],f[2]),(f[2],f[0])):
   directed[(aa,bb)]=directed.get((aa,bb),0)+1;ee=tuple(sorted((aa,bb)));undirected[ee]=undirected.get(ee,0)+1
 expected={tuple(sorted((12+i,12+(i+1)%9))) for i in range(9)};actual=set()
 for ee,nn in undirected.items():
  assert nn in (1,2)
  if nn==1:actual.add(ee)
  else:
   aa,bb=ee;assert directed.get((aa,bb),0)==1 and directed.get((bb,aa),0)==1
 assert actual==expected and len(undirected)==51 and 21-len(undirected)+31==1
 amin=None
 for f in core:
  a=orient(C[f[0]],C[f[1]],C[f[2]]);assert a.l>0;amin=a.lo() if amin is None or a.lo()<amin else amin
 for i in range(9):
  a=orient(Qp[(i-1)%9],Qp[i],C[i]);b=orient(Qp[i],C[(i+1)%9],C[i]);assert a.l>0 and b.l>0
  amin=min(amin,a.lo(),b.lo())
 # boundary convexity/order
 bmin=None
 for i in range(9):
  a=orient(Qp[i],Qp[(i+1)%9],Qp[(i+2)%9]);assert a.l>0;bmin=a.lo() if bmin is None or a.lo()<bmin else bmin
 # circumradii: six active generic/axial faces are exact by defining equations.
 active={(1,2,11),(4,5,9),(7,8,10),(0,10,11),(3,9,11),(6,9,10)}
 inactive=[];active_err=[]
 for fs in [(0,1,11),(0,8,10),(0,10,11),(1,2,11),(2,3,11),(3,4,9),(3,9,11),(4,5,9),(5,6,9),(6,7,10),(6,9,10),(7,8,10),(9,10,11)]:
  # orient locally for circumcenter
  a,b,c0=[C[j] for j in fs]
  if orient(a,b,c0).u<0:b,c0=c0,b
  R2=circum2(a,b,c0)
  if fs in active:
   # enclosure must overlap exact r^2; equality follows symbolically from root equations.
   rr=r*r;assert R2.l<=rr.u and rr.l<=R2.u;active_err.append(float(max(abs(R2.lo()-rr.hi()),abs(R2.hi()-rr.lo()))))
  else:
   assert R2.hi()<(r*r).lo(),(fs,float(R2.hi()-(r*r).lo()));inactive.append((fs,float((r*r).lo()-R2.hi())))
 report={'verified':True,'root_radius':str(rho),'fixed_bits':BITS,'vertices':21,'edges':51,'faces':31,'boundary_edges':9,'combinatorial_disk_verified':True,'min_oriented_area':float(amin),'min_boundary_convexity':float(bmin),'min_inactive_r2_margin':min(x[1] for x in inactive),'inactive_faces':[{'face':list(f),'r2_margin':q} for f,q in inactive],'max_active_interval_error':max(active_err),'gap_half_sine_margin_generic':float(sinDl-s.hi()),'gap_half_sine_margin_axial':float(sinDl-M.hi()),'seconds':time.time()-t}
 json.dump(report,open(ROOT/'cover12_upper_exact_verified.json','w'),indent=2);print(json.dumps(report,indent=2));return report
if __name__=='__main__':verify()
