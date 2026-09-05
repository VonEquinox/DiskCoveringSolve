"""Bind independently audited enumeration, exact metric cuts and residual graphs."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json,struct,hashlib,time
from exact_arcs import sincos_turn,SCALE
from exact_force16 import signature,A,C,C6,D,RU
import system16 as SY
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration';AUDIT=R.parent/'reports';RUNTIME=R.parent/'_runtime';SCALE_LP=10**12

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
@lru_cache(None)
def labelled(B,I):
 assert B>=2 and I>=0
 if B==2:return int(I==0)
 total=sum(comb(I,j)*labelled(k,j)*labelled(B-k+1,I-j) for k in range(2,B) for j in range(I+1))
 if I:total+=I*(labelled(B+1,I-1)-sum(comb(I-1,j)*labelled(3,j)*labelled(B,I-1-j) for j in range(I)))
 return total

def verify_lp(B,key,nums):
 masks=[1<<i for i in range(B)]+list(key)
 assert len(nums)==len(masks) and all(type(v)is int and v>=0 for v in nums)
 assert all(sum(v for v,m in zip(nums,masks) if m>>i&1)>=SCALE_LP for i in range(B))
 assert all(type(m) is int and 0<m<(2<<B) and m&((1<<B)-1) for m in key)
 assert tuple(sorted(set(key)))==tuple(key)
 cost=A*sum(nums[:B])+sum(v*(C6 if m>>B else C) for v,m in zip(nums[B:],key));assert cost<D*SCALE_LP
 return D*SCALE_LP-cost

def check_constants():
 d=json.load(open(R/'root16.json'));t=F(int(d['xnum'][SY.TID]),int(d['Qx']));rho=F(1,10**80)
 assert 0<t-rho<t+rho<RU**2<F(1,4)
 _,s1=sincos_turn(F(A,2*D));_,s2=sincos_turn(F(C,2*D));_,s3=sincos_turn(F(C6,2*D));_,s10=sincos_turn(F(1,20))
 assert RU<F(s1[0],SCALE) and 2*RU<F(s2[0],SCALE) and 3*RU<F(s3[0],SCALE) and RU<F(s10[0],SCALE)
 assert 0<2*A<C<C6<D//2 and 3*A<C6
 return {'radius_upper':str(RU),'gap_cap_turns':str(F(A,D)),'four_rod_cap_turns':str(F(C,D)),'six_rod_cap_turns':str(F(C6,D)),'minimum_boundary_cells':11,'sine_margins':[str(F(s1[0],SCALE)-RU),str(F(s2[0],SCALE)-2*RU),str(F(s3[0],SCALE)-3*RU),str(F(s10[0],SCALE)-RU)]}

def verify():
 st=time.time();const=check_constants();cert=json.load(open(ENUM/'metric16_certificate.json'));rows=json.load(open(ENUM/'metric16_residuals.json'));summary={};tot=0;mass={};cursor=0;minlp=None
 assert set(cert)=={f'B{B}_I{16-B}' for B in range(15,10,-1)}
 for B in range(15,10,-1):
  fam=f'B{B}_I{16-B}';c=cert[fam];maskfile=RUNTIME/(fam+'.masks');assert sha(maskfile)==c['source_masks_sha256']
  ar=json.load(open(AUDIT/f'audit16_B{B}.json'));assert ar['verified'] is True and ar['B']==B and ar['I']==16-B
  expected=labelled(B,16-B);assert int(ar['orbit_mass'])==expected;mass[fam]=expected
  cache={tuple(map(int,k.split(','))) if k else ():v for k,v in c['certified_signatures'].items()};survive={tuple(map(int,k.split(','))) if k else () for k in c['surviving_signatures']};assert not(set(cache)&survive)
  for key,v in cache.items():
   mar=verify_lp(B,key,v);minlp=mar if minlp is None else min(minlp,mar)
  counts={'total':ar['representatives'],'direct':ar['direct'],'farkas':ar['farkas'],'residual':ar['residual']};assert counts==c['counts'] and counts['total']==sum(counts[k] for k in ['direct','farkas','residual'])
  inds=[];seen_survive=set();p=AUDIT/f'audit16_B{B}.residuals'
  with p.open('rb') as f:
   head=f.read(20);assert len(head)==20
   magic,b,i,nf,nr=struct.unpack('<5I',head);assert (magic,b,i,nf,nr)==(0x52534544,B,16-B,30-B,counts['residual'])
   for _ in range(nr):
    raw=f.read(4+2*nf);assert len(raw)==4+2*nf;idx=struct.unpack_from('<I',raw)[0];masks=struct.unpack_from('<'+str(nf)+'H',raw,4)
    assert not inds or idx>inds[-1];assert idx<counts['total'];inds.append(idx)
    faces=[[j for j in range(16) if m>>j&1] for m in masks];assert all(len(ff)==3 for ff in faces)
    sig=signature(B,faces);assert sig is not None and sig in survive and sig not in cache;seen_survive.add(sig)
    row={'B':B,'I':16-B,'idx':idx,'faces':faces,'signature':list(sig)}
    assert cursor<len(rows) and row==rows[cursor];cursor+=1
   assert not f.read(1)
  assert inds==c['residual_indices'] and seen_survive==survive
  summary[fam]=counts;tot+=counts['total'];print('PARTITION FAMILY VERIFIED',fam,counts,flush=True)
 assert cursor==len(rows) and tot==sum(v['total'] for v in summary.values())
 assert len(rows)>0 and len(summary)==5
 out={'verified':True,'constants':const,'families':summary,'total':tot,'residual':len(rows),'minimum_Farkas_integer_margin':minlp,'labelled_totals_from_recurrence':mass,'root_sha256':sha(R/'root16.json'),'residual_sha256':sha(ENUM/'metric16_residuals.json'),'seconds':time.time()-st}
 (R/'partition16_verified.json').write_text(json.dumps(out,indent=2));print('ALL PARTITION VERIFIED',tot,'residual',len(rows),'seconds',time.time()-st,flush=True);return out
if __name__=='__main__':verify()
