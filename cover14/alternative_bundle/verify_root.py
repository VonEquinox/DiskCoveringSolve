#!/usr/bin/env python3
"""Exact contraction verification of the full 116-variable KKT system."""
from pathlib import Path
from fractions import Fraction as F
import json,time,hashlib
from framework import *
if not __debug__:raise RuntimeError('Do not disable assertions')
BASE=Path(__file__).resolve().parent

def load():
 c=json.loads((BASE/'root_certificate.json').read_text());assert c['schema']=='r14-kkt-rational-root-v1' and c['N']==N
 qx,qy=int(c['Qx']),int(c['Qy']);rho=F(int(c['rho_num']),int(c['rho_den']));assert qx>0 and qy>0 and rho>0
 x=[F(int(v),qx) for v in c['xnum']];Y=[[F(int(v),qy) for v in row] for row in c['Ynum']]
 assert len(x)==N and len(Y)==N and all(len(row)==N for row in Y)
 return c,x,Y,rho

def verify():
 st=time.time();c,x,Y,rho=load();fx,J=values_jac(x);Dr=j_lipschitz_rows();cols=[[(k,J[k][j]) for k in range(N) if J[k][j]] for j in range(N)]
 maxratio=F(0);maxcontr=F(0);large_rho=F(1,10**50);large_ratio=F(0)
 for i in range(N):
  row=Y[i];yf=sum((row[k]*fx[k] for k in range(N)),F(0))
  erow=sum((abs(F(i==j)-sum((row[k]*v for k,v in cols[j]),F(0))) for j in range(N)),F(0))
  q=erow+rho*sum((abs(row[k])*Dr[k] for k in range(N)),F(0));rat=abs(yf)/rho+q
  maxratio=max(maxratio,rat);maxcontr=max(maxcontr,q)
  large_ratio=max(large_ratio,abs(yf)/large_rho+erow+large_rho*sum((abs(row[k])*Dr[k] for k in range(N)),F(0)))
 assert maxratio<1 and maxcontr<1 and large_ratio<1
 assert min(x[WSTART:MUSTART])-rho>0
 assert max(x[MUSTART:])+rho<0
 tlo=x[TID]-rho;thi=x[TID]+rho;assert tlo>0
 # Decimal endpoints are exact rational inputs, verified by squaring.
 rlo=F('0.331732034276234122557815548802409538270942754829017478062828')
 rhi=F('0.331732034276234122557815548802409538270942754829017478062829')
 Rlo=F('3.014481257988179142001870494819761487186097134743704272300777')
 Rhi=F('3.014481257988179142001870494819761487186097134743704272300778')
 assert rlo*rlo<tlo<thi<rhi*rhi
 assert Rlo*Rlo*thi<1<Rhi*Rhi*tlo
 out={'verified':True,'N':N,'large_uniqueness_rho':str(large_rho),'large_box_ratio':str(large_ratio),'t_interval':[str(tlo),str(thi)],'r_strict_interval':[str(rlo),str(rhi)],'R_strict_interval':[str(Rlo),str(Rhi)],'contraction':str(maxcontr),'contraction_float':float(maxcontr),'max_ratio':str(maxratio),'max_ratio_float':float(maxratio),'positive_rod_min':str(min(x[WSTART:MUSTART])-rho),'negative_circle_max':str(max(x[MUSTART:])+rho),'root_sha256':hashlib.sha256((BASE/'root_certificate.json').read_bytes()).hexdigest(),'seconds':time.time()-st}
 (BASE/'root_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':verify()
