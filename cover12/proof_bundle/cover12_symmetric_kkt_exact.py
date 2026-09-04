#!/usr/bin/env python3
"""Exact rational Krawczyk certificate for the D3-symmetric n=12 candidate.

The 14 variables are
  r,u,v,w,c,s,z,P,Q,S,A,B,C,D,
where c=cos(h), s=sin(h), z=sqrt(3), and P,...,D are the seven
D3-orbit stress weights.  The equations combine the five geometric contact
relations, z^2=3, free-node equilibrium, anchor tangential stationarity, and
normalization of the stress.

Generation uses high precision numerics only to discover a rational center and
an approximate inverse. Verification is pure Fraction interval arithmetic.
"""
from __future__ import annotations
from fractions import Fraction as F
import json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
N=14
NAMES=['r','u','v','w','c','s','z','P','Q','S','A','B','C','D']

# Polynomial represented as list (coefficient, exponent tuple).
def addterm(poly,coef,**kw):
 e=[0]*N
 for k,v in kw.items():e[NAMES.index(k)]=v
 poly.append((F(coef),tuple(e)))
def polys():
 r,u,v,w,c,s,z,P,Q,S,A,B,C,D=range(N);out=[]
 def p(*terms):out.append([(F(a),tuple(e)) for a,e in terms])
 def E(**kw):
  q=[0]*N
  for k,vv in kw.items():q[NAMES.index(k)]=vv
  return tuple(q)
 # c^2+s^2-1
 p((1,E(c=2)),(1,E(s=2)),(-1,E()))
 # u^2-2uc+1-r^2
 p((1,E(u=2)),(-2,E(u=1,c=1)),(1,E()),(-1,E(r=2)))
 # v^2-v(c^2-s^2+2zcs)+1-r^2
 p((1,E(v=2)),(-1,E(v=1,c=2)),(1,E(v=1,s=2)),(-2,E(v=1,z=1,c=1,s=1)),(1,E()),(-1,E(r=2)))
 # w-u^2+r^2+r
 p((1,E(w=1)),(-1,E(u=2)),(1,E(r=2)),(1,E(r=1)))
 # v^2-(2r+w)v+w(w+r)
 p((1,E(v=2)),(-2,E(r=1,v=1)),(-1,E(w=1,v=1)),(1,E(w=2)),(1,E(w=1,r=1)))
 # z^2-3
 p((1,E(z=2)),(-3,E()))
 # 2A(uc-(w+r))-Br
 p((2,E(A=1,u=1,c=1)),(-2,E(A=1,w=1)),(-2,E(A=1,r=1)),(-1,E(B=1,r=1)))
 # Cr-2D(v-r-w/2) = Cr-2Dv+2Dr+Dw
 p((1,E(C=1,r=1)),(-2,E(D=1,v=1)),(2,E(D=1,r=1)),(1,E(D=1,w=1)))
 # 2P(L-v)-Cr, L=(c^2-s^2)/2+zcs
 p((1,E(P=1,c=2)),(-1,E(P=1,s=2)),(2,E(P=1,z=1,c=1,s=1)),(-2,E(P=1,v=1)),(-1,E(C=1,r=1)))
 # Q-S-A(w+r)
 p((1,E(Q=1)),(-1,E(S=1)),(-1,E(A=1,w=1)),(-1,E(A=1,r=1)))
 # (S+Q)(c-u)+A((w+r)c-u)
 p((1,E(S=1,c=1)),(-1,E(S=1,u=1)),(1,E(Q=1,c=1)),(-1,E(Q=1,u=1)),(1,E(A=1,w=1,c=1)),(1,E(A=1,r=1,c=1)),(-1,E(A=1,u=1)))
 # Br-D(2w-v+r)
 p((1,E(B=1,r=1)),(-2,E(D=1,w=1)),(1,E(D=1,v=1)),(-1,E(D=1,r=1)))
 # Qus-PvM; M=z(c^2-s^2)/2-cs
 p((1,E(Q=1,u=1,s=1)),(-F(1,2),E(P=1,v=1,z=1,c=2)),(F(1,2),E(P=1,v=1,z=1,s=2)),(1,E(P=1,v=1,c=1,s=1)))
 # normalization
 p((6,E(P=1)),(6,E(Q=1)),(6,E(S=1)),(6,E(A=1)),(6,E(D=1)),(3,E(B=1)),(3,E(C=1)),(-1,E()))
 assert len(out)==N
 return out
POLYS=polys()

def deriv(poly,j):
 out=[]
 for a,e in poly:
  if e[j]:
   ee=list(e);k=ee[j];ee[j]-=1;out.append((a*k,tuple(ee)))
 return out
JAC=[[deriv(POLYS[i],j) for j in range(N)] for i in range(N)]

class Iv:
 __slots__=('l','u')
 def __init__(self,l,u=None):self.l=F(l);self.u=F(l if u is None else u);assert self.l<=self.u
 def __add__(self,o):o=iv(o);return Iv(self.l+o.l,self.u+o.u)
 __radd__=__add__
 def __neg__(self):return Iv(-self.u,-self.l)
 def __sub__(self,o):return self+(-iv(o))
 def __rsub__(self,o):return iv(o)-self
 def __mul__(self,o):
  o=iv(o);q=[self.l*o.l,self.l*o.u,self.u*o.l,self.u*o.u];return Iv(min(q),max(q))
 __rmul__=__mul__
 def __pow__(self,n):
  assert n>=0
  if n==0:return Iv(1)
  if n==1:return self
  if n%2==0:
   if self.l<=0<=self.u:lo=F(0)
   else:lo=min(self.l**n,self.u**n)
   return Iv(lo,max(self.l**n,self.u**n))
  return Iv(self.l**n,self.u**n)
 def absup(self):return max(abs(self.l),abs(self.u))
 def width(self):return self.u-self.l
 def __repr__(self):return f'[{float(self.l):.4g},{float(self.u):.4g}]'
def iv(x):return x if isinstance(x,Iv) else Iv(x)

def peval(poly,x):
 s=Iv(0)
 for a,e in poly:
  q=Iv(a)
  for j,k in enumerate(e):
   if k:q=q*(x[j]**k)
  s=s+q
 return s

def point_eval(poly,x):
 s=F(0)
 for a,e in poly:
  q=a
  for j,k in enumerate(e):
   if k:q*=x[j]**k
  s+=q
 return s

def matmul(A,B):
 return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Iv(0)) for j in range(len(B[0]))] for i in range(len(A))]

def generate(path):
 import sympy as sp, mpmath as mp
 sy=sp.symbols(' '.join(NAMES));r,u,v,w,c,s,z,P,Q,S,A,B,C,D=sy
 L=(c*c-s*s)/2+z*c*s;M=z*(c*c-s*s)/2-c*s
 eq=[c*c+s*s-1,u*u-2*u*c+1-r*r,v*v-2*v*L+1-r*r,w-u*u+r*r+r,v*v-(2*r+w)*v+w*(w+r),z*z-3,
 2*A*(u*c-(w+r))-B*r,C*r-2*D*(v-r-w/2),2*P*(L-v)-C*r,Q-S-A*(w+r),(S+Q)*(c-u)+A*((w+r)*c-u),B*r-D*(2*w-v+r),Q*u*s-P*v*M,6*(P+Q+S+A+D)+3*(B+C)-1]
 base=json.load(open(ROOT/'cover12_candidate_symbolic.json'));h=sp.Float(base['h'],150);st=base['stress'];guess=[base['r'],base['u'],base['v'],base['w'],sp.cos(h),sp.sin(h),sp.N(sp.sqrt(3),150),st[0],st[1],st[2],st[18],st[20],st[27],st[28]]
 sol=sp.nsolve(eq,sy,guess,tol=sp.Float('1e-135'),maxsteps=100,prec=160)
 Qx=10**120;xn=[]
 for q in sol:
  a=F(str(sp.N(q,150)));zz=a*Qx;xn.append((2*zz.numerator+zz.denominator)//(2*zz.denominator) if zz>=0 else -((2*(-zz.numerator)+zz.denominator)//(2*zz.denominator)))
 xc=[mp.mpf(n)/Qx for n in xn];mp.mp.dps=180
 # Jacobian via exact polynomials evaluated in mp
 def mpeval(poly):
  ss=mp.mpf('0')
  for a,e in poly:
   q=mp.mpf(a.numerator)/a.denominator
   for j,k in enumerate(e):q*=xc[j]**k
   ss+=q
  return ss
 Jm=mp.matrix([[mpeval(JAC[i][j]) for j in range(N)] for i in range(N)]);Ym=Jm**-1;Qy=10**100
 yn=[[int(mp.nint(Ym[i,j]*Qy)) for j in range(N)] for i in range(N)]
 cert={'version':1,'names':NAMES,'Qx':str(Qx),'xnum':[str(n) for n in xn],'Qy':str(Qy),'Ynum':[[str(n) for n in row] for row in yn],'rho_num':'1','rho_den':str(10**80)}
 json.dump(cert,open(path,'w'),separators=(',',':'));print('saved',path);verify(path)

def verify(path):
 t0=time.time();d=json.load(open(path));assert d['version']==1 and d['names']==NAMES
 Qx=int(d['Qx']);x=[F(int(n),Qx) for n in d['xnum']];Qy=int(d['Qy']);Y=[[F(int(n),Qy) for n in row] for row in d['Ynum']];rho=F(int(d['rho_num']),int(d['rho_den']))
 assert len(x)==N and len(Y)==N and all(len(row)==N for row in Y)
 X=[Iv(q-rho,q+rho) for q in x];f=[point_eval(p,x) for p in POLYS];J=[[peval(JAC[i][j],X) for j in range(N)] for i in range(N)]
 # Krawczyk relative interval -Yf+(I-YJ)[-rho,rho]
 YI=[[Iv(q) for q in row] for row in Y];YJ=matmul(YI,J);rat=[];K=[]
 # A separate contraction check makes the Krawczyk uniqueness hypothesis
 # explicit and also proves that both Y and every Jacobian in the box are nonsingular.
 M=[[Iv(1 if i==j else 0)-YJ[i][j] for j in range(N)] for i in range(N)]
 contraction=max(sum((q.absup() for q in row),F(0)) for row in M)
 assert contraction<1
 for i in range(N):
  cf=-sum((Y[i][k]*f[k] for k in range(N)),F(0));z=Iv(cf)
  for j in range(N):
   m=M[i][j];z=z+m*Iv(-rho,rho)
  assert -rho<z.l and z.u<rho,(i,float(z.l/rho),float(z.u/rho))
  rat.append(max(abs(z.l),abs(z.u))/rho);K.append([z.l,z.u])
 # positivity and ordering/branch inequalities
 assert all(q-rho>0 for q in x)
 # h branch represented by c,s: c,s>0 and 2h<pi/3 follows later via algebraic inequality.
 report={'verified':True,'max_krawczyk_ratio':float(max(rat)),'contraction_inf_norm':float(contraction),'rho':str(rho),'min_variable_lower':float(min(q-rho for q in x)),'r_interval':[str(x[0]-rho),str(x[0]+rho)],'t_interval':[str((x[0]-rho)**2),str((x[0]+rho)**2)],'seconds':time.time()-t0,'center_decimal':{n:format(float(q),'.17g') for n,q in zip(NAMES,x)}}
 json.dump(report,open(ROOT/'cover12_symmetric_kkt_verified.json','w'),indent=2);print(json.dumps(report,indent=2))
 return report
if __name__=='__main__':
 mode=sys.argv[1] if len(sys.argv)>1 else 'verify';path=sys.argv[2] if len(sys.argv)>2 else str(ROOT/'cover12_symmetric_kkt_cert.json')
 generate(path) if mode=='generate' else verify(path)
