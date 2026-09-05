#!/usr/bin/env python3
"""Exact polynomial/graph/enum/root interfaces and decimal isolation."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from functools import lru_cache
from math import comb,factorial
import json
if not __debug__:
 raise RuntimeError('Exact checks require Python without -O/PYTHONOPTIMIZE')
ROOT=Path(__file__).resolve().parent;ENUM=ROOT.parent/'enumeration'

@lru_cache(None)
def triangulation_count(b,i):
 """Root-edge peeling recurrence; boundary labels fixed, i interior labels distinct."""
 if b==2:return int(i==0)
 assert b>=3 and i>=0
 s=sum(comb(i,j)*triangulation_count(k,j)*triangulation_count(b-k+1,i-j) for k in range(2,b) for j in range(i+1))
 if i:
  forbidden=sum(comb(i-1,j)*triangulation_count(3,j)*triangulation_count(b,i-1-j) for j in range(i))
  available=triangulation_count(b+1,i-1)-forbidden
  assert available>=0;s+=i*available
 return s

def verify():
 import verify_root_repaired as K
 L=json.loads((ROOT/'cover13_core_layout.json').read_text());rd=json.loads((ROOT/'cover13_kkt_root_119d.json').read_text());cert=json.loads((ROOT/'cover13_krawczyk_cert.json').read_text())
 assert L['layout']=='qcp' and L['constraint_order']=='rn' and L['lambda_sign']==1 and L['objective_sign']==1
 assert L['constraint_labels']==[['rod',i] for i in range(47)]+[['qnorm',i] for i in range(1,10)]
 cons,ee=K.build(L)
 def node(v):
  t,i=v;return i if t=='q' else 10+i if t=='c' else 23+i if t=='p' else None
 ed=[[node(u),node(v)] for u,v in ee]
 assert ed==rd['edges']
 # Every polynomial is exactly the squared-distance/unit-circle equation in this canonical graph.
 assert ed[:20]==[[u,v] for i in range(10) for u,v in [(i,10+i),(i,10+(i+1)%10)]]
 assert ed[20:]==[[23+k,10+v] for k,fa in enumerate(rd['active_faces']) for v in fa]
 assert len(rd['active_faces'])==9 and len(rd['all_faces'])==14
 assert set(map(tuple,rd['active_faces']))<=set(map(tuple,rd['all_faces']))
 cm=json.loads((ENUM/'candidate_orbit_map.json').read_text());rr=json.loads((ENUM/'metric_residuals_new.json').read_text());k=cm['residual_index'];row=rr[k]
 assert (row['B'],row['I'],row['idx'])==(10,3,cm['enumeration_idx'])
 refl,shift,p=cm['map'];assert refl in (0,1) and 0<=shift<10 and sorted(p)==list(range(13))
 assert all(p[i]==((shift-i) if refl else (shift+i))%10 for i in range(10))
 assert sorted(p[10:])==list(range(10,13))
 assert {tuple(sorted(p[v] for v in fa)) for fa in rd['all_faces']}==set(map(tuple,row['faces']))
 counts={}
 for B in range(9,13):
  I=13-B;v=triangulation_count(B,I)
  num=2*factorial(2*B-3)*factorial(4*I+2*B-5);den=factorial(B-1)*factorial(B-3)*factorial(3*I+2*B-3)
  assert num%den==0 and v==num//den
  counts[f'B{B}_I{I}']=str(v)
 t=F(int(cert['xnum'][62]),int(cert['Qx']));rho=F(int(cert['rho_num']),int(cert['rho_den']))
 assert rho>0 and t>rho
 # Decimal calculations propose endpoints only; the following integer comparisons accept them.
 with localcontext() as c:
  c.prec=100;dv=(Decimal(t.numerator)/Decimal(t.denominator)).sqrt();scale=10**60;z=int(dv*scale)
  rv=1/dv;zz=int(rv*scale)
 rlo=F(z,scale);rhi=F(z+1,scale);Rlo=F(zz,scale);Rhi=F(zz+1,scale)
 assert 0<rlo and rlo*rlo<t-rho<t+rho<rhi*rhi
 assert 0<Rlo and Rlo*Rlo*(t+rho)<1 and Rhi*Rhi*(t-rho)>1
 def dec(n):return str(n//scale)+'.'+str(n%scale).zfill(60)
 out={'verified':True,'dimension':119,'active_rods':47,'all_center_faces':14,'candidate_residual_index':k,'candidate_enumeration_index':row['idx'],'candidate_permutation':p,'labeled_counts_by_recurrence':counts,'r_strict_interval':[dec(z),dec(z+1)],'R_strict_interval':[dec(zz),dec(zz+1)],'t_root_interval':[str(t-rho),str(t+rho)]}
 (ROOT/'interfaces_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':verify()
