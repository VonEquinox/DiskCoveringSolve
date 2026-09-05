"""Rebuild the integer KKT system and accept two nested rational contraction boxes."""
from fractions import Fraction as F
import json,time,hashlib
from math import isqrt
from pathlib import Path
import system16 as S
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent

def verify():
 st=time.time();d,z,Y,rho=S.load();assert d['N']==S.D and rho==F(1,10**60);assert len(z)==S.D and len(Y)==S.D and all(len(r)==S.D for r in Y)
 assert d['edges']==[list(e) for e in S.EDGES] and d['active_faces']==[list(e) for e in S.ACTIVE]
 f,J=S.evaluate(z,F(0));n=S.D
 YF=[sum((a*b for a,b in zip(row,f)),F(0)) for row in Y]
 jc=[[(k,J[k][j]) for k in range(n) if J[k][j]] for j in range(n)]
 Erow=[]
 for i in range(n):
  Erow.append(sum(abs(F(int(i==j))-sum((Y[i][k]*v for k,v in jc[j]),F(0))) for j in range(n)))
 Dr=[F(0)]*n
 for i,(H,b,c) in enumerate(S.CONS):
  Dr[i]=sum(abs(v) for v in H.values())
  for (j,k),v in H.items():Dr[S.C+j]+=2*abs(v)
 der=[sum((abs(a)*v for a,v in zip(row,Dr)),F(0)) for row in Y]
 reports=[]
 for r in [F(1,10**60),F(1,10**80)]:
  ratios=[abs(a)/r+e+r*h for a,e,h in zip(YF,Erow,der)];qs=[e+r*h for e,h in zip(Erow,der)]
  assert max(ratios)<1 and max(qs)<1
  reports.append({'rho':str(r),'inclusion_ratio':str(max(ratios)),'contraction':str(max(qs))})
 assert all(v-F(1,10**80)>0 for v in z[S.P:S.P+S.M]);assert all(-F(11,1000)<v-F(1,10**80) and v+F(1,10**80)<0 for v in z[S.P+S.M:])
 digits=60;scale=10**digits;tlo=z[S.TID]-F(1,10**80);thi=z[S.TID]+F(1,10**80);tm=z[S.TID]
 k=isqrt(tm.numerator*scale*scale//tm.denominator);kr=isqrt(tm.denominator*scale*scale//tm.numerator)
 rlo=F(k,scale);rhi=F(k+1,scale);Rlo=F(kr,scale);Rhi=F(kr+1,scale)
 assert rlo*rlo<tlo<thi<rhi*rhi
 assert Rlo*Rlo*thi<1<Rhi*Rhi*tlo
 def dec(j):return str(j//scale)+'.'+str(j%scale).zfill(digits)
 decimals={'r_lower':dec(k),'r_upper':dec(k+1),'R_lower':dec(kr),'R_upper':dec(kr+1)}
 out={'verified':True,'dimension':n,'boxes':reports,'minimum_weight':str(min(z[S.P:S.P+S.M])-F(1,10**80)),'r_interval':[str(rlo),str(rhi)],'R_interval':[str(Rlo),str(Rhi)],'decimal_intervals':decimals,'root_sha256':hashlib.sha256((R/'root16.json').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (R/'root16_verified.json').write_text(json.dumps(out,indent=2));print('ROOT VERIFIED',n,'ratio',float(F(reports[1]['inclusion_ratio'])),'r',float(rlo),'sec',time.time()-st,flush=True)
 return out
if __name__=='__main__':verify()
