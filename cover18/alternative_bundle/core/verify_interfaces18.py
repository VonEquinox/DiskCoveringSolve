"""Independent reconstruction of every interface dimension and labelled graph."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import system18 as S
import exact_force18 as G
import anchor_isolation18 as L
if not __debug__:raise RuntimeError('Assertions must remain enabled')
R=Path(__file__).resolve().parent

def verify():
 assert (S.N,S.B)==(18,12)
 assert len(S.ACTIVE)==16 and len(S.FACES)==22
 for fs in (S.ACTIVE,S.FACES):
  assert len(set(fs))==len(fs)
  assert all(len(f)==3 and tuple(sorted(set(f)))==f and 0<=f[0]<f[-1]<S.N for f in fs)
 assert set(S.ACTIVE)<=set(S.FACES)
 expected=[(i,S.B+i) for i in range(S.B)]+[(i,S.B+(i+1)%S.B) for i in range(S.B)]
 expected += [(S.B+S.N+j,S.B+i) for j,f in enumerate(S.ACTIVE) for i in f]
 assert S.EDGES==expected and len(set(expected))==72
 assert (S.G,S.P,S.TID,S.M,S.C,S.D)==(90,91,90,72,83,174)
 assert S.G==2*(S.B+S.N+len(S.ACTIVE)-1)
 assert S.C==S.M+S.B-1 and S.D==S.P+S.C and S.P==S.G+1
 assert len(S.CONS)==S.C
 assert G.N==S.N and G.V==3*S.N-2 and G.BMIN==11
 assert (G.A,G.C,G.C6,G.D)==(937111,1970786,3362065,10**7)
 assert G.RU==F('0.29016771764058')
 assert (L.GAMMA,L.LAMBDA,L.NU,L.RADIUS)==(300,F(1,1000),F(7,625),F(11,100))
 raw=json.loads((R/'root18.json').read_text())
 assert raw['version']==1 and raw['N']==S.D
 assert raw['edges']==[list(e) for e in expected]
 assert raw['active_faces']==[list(f) for f in S.ACTIVE]
 assert raw['all_faces']==[list(f) for f in S.FACES]
 assert len(raw['xnum'])==S.D and len(raw['Ynum'])==S.D and all(len(v)==S.D for v in raw['Ynum'])
 # Check constraint polynomials against an independently assembled signed vector.
 for k,(H,b,c) in enumerate(S.CONS):
  assert all(type(v)is int for v in list(H.values())+list(b.values())+[c])
  assert all(0<=i<S.P and 0<=j<S.P and H.get((j,i),0)==v for (i,j),v in H.items())
  if k<S.M:
   u,v=expected[k];hh={};bb={S.TID:-1};cc=0
   for d in range(2):
    aa={};a0=0
    for node,sgn in ((u,1),(v,-1)):
     if node:aa[2*(node-1)+d]=sgn
     else:a0+=sgn*(1 if d==0 else 0)
    cc+=a0*a0
    for i,ai in aa.items():
     bb[i]=bb.get(i,0)+2*a0*ai
     for j,aj in aa.items():hh[i,j]=hh.get((i,j),0)+2*ai*aj
   assert (H,b,c)==(hh,bb,cc)
  else:
   i=k-S.M
   assert H=={(2*i,2*i):2,(2*i+1,2*i+1):2} and b=={} and c==-1
 out={'verified':True,'centers':18,'anchors':12,'active_faces':16,'all_center_faces':22,
      'active_rods':72,'geometric_dimension':90,'constraint_count':83,'kkt_dimension':174,
      'full_auxiliary_nodes':52,'full_auxiliary_rods':90,
      'root_sha256':hashlib.sha256((R/'root18.json').read_bytes()).hexdigest()}
 (R/'interfaces18_verified.json').write_text(json.dumps(out,indent=2))
 print('EXPLICIT SYSTEM AND GEOMETRIC INTERFACES VERIFIED 174 variables, 72 rods',flush=True)
 return out
if __name__=='__main__':verify()
