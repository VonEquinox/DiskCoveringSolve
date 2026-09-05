"""Exact, nonnumerical cross-stage linkage and analytic pruning premises."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
import system18 as S
import exact_force18 as G
import anchor_isolation18 as L
from exact_arcs import sincos_turn,SCALE,QA,QB
if not __debug__:raise RuntimeError('Run without -O/PYTHONOPTIMIZE')
R=Path(__file__).resolve().parent

def verify():
 (R/'interfaces18_verified.json').unlink(missing_ok=True)
 st=time.time()
 assert (S.N,S.B,S.G,S.P,S.M,S.C,S.D,S.TID)==(18,12,90,91,72,83,174,90)
 assert len(S.FACES)==22 and len(set(S.FACES))==22
 assert len(S.ACTIVE)==16 and set(S.ACTIVE)<=set(S.FACES)
 expected=[(i,S.B+i) for i in range(S.B)]+[(i,S.B+(i+1)%S.B) for i in range(S.B)]
 expected +=[(S.B+S.N+j,S.B+v) for j,f in enumerate(S.ACTIVE) for v in f]
 assert expected==S.EDGES and len(S.CONS)==S.C
 assert all(tuple(sorted(f))==f and len(set(f))==3 and min(f)>=0 and max(f)<S.N for f in S.FACES)
 assert (G.N,G.V,G.BMIN)==(18,52,11)
 assert (G.A,G.C,G.C6,G.D)==(9372,19708,33621,100000)
 assert 0<G.A<G.C<G.C6<G.D//2 and G.RU==F('0.29016771764058')
 assert QA==G.QL==10**15 and QB==10**18 and G.DF%(G.QL*QA)==0
 _,z,_,_=S.load();eps=F(1,10**80)
 assert 0<z[S.TID]-eps<z[S.TID]+eps<G.RU**2
 margins={}
 for key,cap,mul in [('two_rod',G.A,1),('four_rod',G.C,2),('six_rod',G.C6,3)]:
  co,si=sincos_turn(F(cap,2*G.D));margin=F(si[0],SCALE)-mul*G.RU;assert margin>0;margins[key]=str(margin)
 _,si=sincos_turn(F(1,20));assert G.RU<F(si[0],SCALE) and G.RU<F(1,2)
 assert L.RADIUS==F(7,60) and L.GAMMA==400 and L.LAMBDA==F(69,50000) and L.NU==F(7,625)
 assert L.LAMBDA-2*L.GAMMA*L.NU**2*L.RADIUS**2==F(793,56250000)>0
 out={'verified':True,'N':S.N,'families':[[b,18-b] for b in range(11,18)],'radius_upper':str(G.RU),'angle_denominator':G.D,'gap_cap':G.A,'four_rod_cap':G.C,'six_rod_cap':G.C6,'sine_margins':margins,'root_sha256':hashlib.sha256((R/'root18.json').read_bytes()).hexdigest(),'anchor_radius':str(L.RADIUS),'seconds':time.time()-st}
 (R/'interfaces18_verified.json').write_text(json.dumps(out,indent=2));print('EXACT INTERFACES VERIFIED',{k:float(F(v)) for k,v in margins.items()},flush=True);return out
if __name__=='__main__':verify()
