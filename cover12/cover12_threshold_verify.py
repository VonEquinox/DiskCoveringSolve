#!/usr/bin/env python3
"""Exact bridge between the algebraic candidate root and all rational thresholds."""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from fixed_trig import sincos_point
from exact_core import PI_L,PI_U
T=F(814971,6250000)       # 0.13039536
DELTA=F(73891,100000)     # boundary arc cap
ALPHA={2:F(73891,100000),3:F(57241,50000),4:F(161399,100000),5:F(22523,10000)}

def verify():
 t0=time.time();d=json.load(open(ROOT/'cover12_symmetric_kkt_cert.json'))
 q=int(d['Qx']);rho=F(int(d['rho_num']),int(d['rho_den']));r0=F(int(d['xnum'][0]),q)
 rl,ru=r0-rho,r0+rho;tl,tu=rl*rl,ru*ru
 assert 0<rl<ru<F(1,2)
 assert tu<T
 sl,su,cl,cu=sincos_point(DELTA/2)
 boundary_margin=sl*sl-T
 assert boundary_margin>0
 assert 8*DELTA<2*PI_L
 assert 4*DELTA<PI_L
 angle={}
 for dd,a in ALPHA.items():
  ss,_,_,_=sincos_point(a/2);mar=4*ss*ss-dd*dd*T;assert mar>0
  angle[str(dd)]=str(mar)
 report={
  'verified':True,'T':str(T),'delta':str(DELTA),
  'r_interval':[str(rl),str(ru)],'t_interval':[str(tl),str(tu)],
  'T_minus_t_upper':str(T-tu),'T_minus_t_upper_float':float(T-tu),
  'half_minus_r_upper':str(F(1,2)-ru),'half_minus_r_upper_float':float(F(1,2)-ru),
  'sin2_delta_half_minus_T':str(boundary_margin),'sin2_delta_half_minus_T_float':float(boundary_margin),
  'two_pi_minus_8delta':str(2*PI_L-8*DELTA),'two_pi_minus_8delta_float':float(2*PI_L-8*DELTA),
  'pi_minus_4delta':str(PI_L-4*DELTA),'pi_minus_4delta_float':float(PI_L-4*DELTA),
  'angle_chord_margins':angle,'seconds':time.time()-t0,
 }
 (ROOT/'cover12_threshold_verified.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':verify()
