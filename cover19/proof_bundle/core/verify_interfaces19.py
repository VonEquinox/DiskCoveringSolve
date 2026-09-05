"""Exact link from 1/sqrt(13) to all necessary-condition constants."""
from pathlib import Path
from fractions import Fraction as F
import json,time
import geometry19 as S
import exact_force19 as G
from exact_arcs import sincos_turn,SCALE,QA,QB
if not __debug__:raise RuntimeError('Assertions must be enabled')
R=Path(__file__).resolve().parent

def verify():
 (R/'interfaces19_verified.json').unlink(missing_ok=True);t=time.time()
 assert (G.N,G.V,G.BMIN,G.D)==(19,55,12,100000)
 assert G.QL==QA and QB%QA==0 and G.DF%(G.QL*QA)==0
 assert S.T==F(1,13) and G.RU>0 and G.RU**2>S.T and G.RU<F(1,2)
 margins=[]
 for k,cap in [(1,G.A),(2,G.C),(3,G.C6)]:
  assert 0<cap<G.D/2
  _,s=sincos_turn(F(cap,2*G.D));mar=F(s[0],SCALE)-k*G.RU;assert mar>0;margins.append(str(mar))
 _,s=sincos_turn(F(1,22));assert F(s[0],SCALE)>G.RU
 assert 2*G.A<G.C<G.C6 and 12*G.A>G.D
 out={'verified':True,'radius_squared':'1/13','radius_upper':str(G.RU),'turn_caps':[str(F(x,G.D)) for x in (G.A,G.C,G.C6)],'sine_margins':margins,'minimum_boundary_cells':12,'families':[[b,19-b] for b in range(12,19)],'seconds':time.time()-t}
 (R/'interfaces19_verified.json').write_text(json.dumps(out,indent=2));print('EXACT GEOMETRY / ENUMERATION INTERFACES VERIFIED',flush=True);return out
if __name__=='__main__':verify()
