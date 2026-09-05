"""Complete exact metric partition, recomputed from independently audited faces."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from functools import lru_cache
from math import comb
import json,gzip,hashlib,time
from exact_arcs import sincos_turn,SCALE
from exact_force14 import signature,A,C,D,RU
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration';SCALE_LP=10**12
@lru_cache(maxsize=None)
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
 cost=A*sum(nums[:B])+C*sum(nums[B:]);assert cost<D*SCALE_LP
 return D*SCALE_LP-cost

def check_constants():
 d=json.load(open(R/'root14.json'));t=F(int(d['xnum'][62]),int(d['Qx']));rho=F(1,10**80)
 assert 0<t-rho<t+rho<RU**2<F(1,4)
 _,s1=sincos_turn(F(A,2*D));_,s2=sincos_turn(F(C,2*D));_,s9=sincos_turn(F(1,18))
 assert RU<F(s1[0],SCALE) and 2*RU<F(s2[0],SCALE) and RU<F(s9[0],SCALE)
 assert 0<2*A<C<D//2
 return {'radius_upper':str(RU),'gap_cap_turns':str(F(A,D)),'four_rod_cap_turns':str(F(C,D)),'minimum_boundary_cells':10,'sine_margins':[str(F(s1[0],SCALE)-RU),str(F(s2[0],SCALE)-2*RU),str(F(s9[0],SCALE)-RU)]}

def verify():
 st=time.time();const=check_constants();cert=json.load(open(ENUM/'metric14_certificate.json'));triples=list(combinations(range(14),3));allres=[];summary={};tot=0;mass={}
 for B in range(10,14):
  fam=f'B{B}_I{14-B}';p=ENUM/(fam+'.txt.gz')
  raw=gzip.decompress(p.read_bytes()) if p.exists() else (ENUM/(fam+'.txt')).read_bytes();c=cert[fam];assert hashlib.sha256(raw).hexdigest()==c['source_txt_sha256']
  lines=raw.decode('ascii').splitlines();b,i,nf,nstates=map(int,lines[0].split());assert (b,i,nf)==(B,14-B,26-B) and len(lines)-1==nstates
  cache={tuple(map(int,k.split(','))) if k else ():v for k,v in c['certified_signatures'].items()};survive={tuple(map(int,k.split(','))) if k else () for k in c['surviving_signatures']}
  assert not(set(cache)&survive)
  for key,v in cache.items():verify_lp(B,key,v)
  counts=Counter({'direct':0,'farkas':0,'residual':0});inds=[];seen_survive=set()
  for idx,line in enumerate(lines[1:]):
   ids=list(map(int,line.split()));assert len(ids)==nf and len(set(ids))==nf and all(0<=a<len(triples) for a in ids)
   faces=[list(triples[a]) for a in ids];key=signature(B,faces)
   if key is None:counts['direct']+=1
   elif key in cache:counts['farkas']+=1
   else:
    assert key in survive;seen_survive.add(key);counts['residual']+=1;inds.append(idx)
    allres.append({'B':B,'I':14-B,'idx':idx,'faces':faces,'signature':list(key)})
  assert dict(counts)==c['counts'] and seen_survive==survive and inds==c['residual_indices'];summary[fam]=dict(counts);tot+=nstates
  audit=json.load(open(R.parent/'reports'/f'audit{B}.json'));assert audit['verified'] and audit['representatives']==nstates and audit['B']==B and audit['I']==14-B
  expected=labelled(B,14-B);assert int(audit['orbit_mass'])==expected;mass[fam]=expected
  print('METRIC VERIFIED',fam,dict(counts),'seconds',time.time()-st,flush=True)
 assert tot==1313024 and len(allres)==7232 and allres==json.load(open(ENUM/'metric14_residuals.json'))
 out={'verified':True,'constants':const,'families':summary,'total':tot,'residual':len(allres),'labelled_totals_from_recurrence':mass,'seconds':time.time()-st}
 (R/'partition14_verified.json').write_text(json.dumps(out,indent=2));print('PARTITION VERIFIED',tot,'residual',len(allres),flush=True);return out
if __name__=='__main__':verify()
