#!/usr/bin/env python3
"""Independent exact angle/Farkas partition replay. Python standard library only."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from functools import lru_cache
from math import comb,factorial
import gzip,json,time,hashlib
from framework import ALL_FACES
from verify_root import load
from exact_arcs import sincos_turn,SCALE
if not __debug__:raise RuntimeError('Assertions must be enabled')
BASE=Path(__file__).resolve().parent;ENUM=BASE/'enumeration'
A=10764;C=23092;D=100000;S=10**12;RU=F(33173203427624,10**14)
FAMILIES=[(10,4),(11,3),(12,2),(13,1)]
def sha(b):return hashlib.sha256(b).hexdigest()
@lru_cache(None)
def pair_rules(B):
 rr=[]
 for i in range(B):
  for j in range(i+1,B):
   k=j-i
   if min(k,B-k)<=2:continue
   l=k*A<D-C;r=(B-k)*A<D-C
   if not(l or r):continue
   arc=((1<<j)-1)^((1<<i)-1)
   mask=None if l and r else arc if l else ((1<<B)-1)^arc
   rr.append((i,j,mask))
 return rr

def signature(B,faces):
 adjacent=[1<<i for i in range(B)]
 for face in faces:
  boundary=[v for v in face if v<B]
  bits=sum(1<<v for v in boundary)
  for v in boundary:adjacent[v]|=bits
 masks=set()
 for i,j,mask in pair_rules(B):
  # A four-rod route exists whenever some incident boundary centers share a face.
  if (adjacent[i]|adjacent[(i+1)%B])&((1<<j)|(1<<((j+1)%B))):
   if mask is None:return None
   masks.add(mask)
 return tuple(sorted(masks))

def verify_lp(B,key,nums):
 masks=[1<<i for i in range(B)]+list(key)
 assert len(nums)==len(masks) and all(type(v)is int and v>=0 for v in nums)
 assert all(sum(n for n,m in zip(nums,masks) if (m>>i)&1)>=S for i in range(B))
 v=D*S-A*sum(nums[:B])-C*sum(nums[B:]);assert v>0;return v

@lru_cache(None)
def triangulation_count(b,i):
 """Root-edge peeling; internal vertices have distinct labels."""
 if b==2:return int(i==0)
 assert b>=3 and i>=0
 total=sum(comb(i,j)*triangulation_count(k,j)*triangulation_count(b-k+1,i-j) for k in range(2,b) for j in range(i+1))
 if i:
  forbidden=sum(comb(i-1,j)*triangulation_count(3,j)*triangulation_count(b,i-1-j) for j in range(i))
  available=triangulation_count(b+1,i-1)-forbidden;assert available>=0;total+=i*available
 return total

def verify():
 st=time.time();cert=json.loads((ENUM/'metric_partition_certificate.json').read_text());given=json.loads((ENUM/'metric_residuals.json').read_text())
 assert set(cert)=={f'B{B}_I{I}' for B,I in FAMILIES}
 c,x,Y,rho=load();assert 0<x[62]-rho<x[62]+rho<RU**2 and 0<RU<F(1,2)
 assert 2*A<C and 2*C<D
 sins={name:F(sincos_turn(v)[1][0],SCALE) for name,v in [('boundary',F(A,2*D)),('four_rod',F(C,2*D)),('nine_boundary',F(1,18))]}
 assert sins['boundary']>RU and sins['four_rod']>2*RU and sins['nine_boundary']>RU
 triples=list(combinations(range(14),3));report={};regenerated=[];minlp=None
 for B,I in FAMILIES:
  name=f'B{B}_I{I}';fc=cert[name];data=(ENUM/(name+'.txt.gz')).read_bytes();raw=gzip.decompress(data)
  assert sha(data)==fc['source_gzip_sha256'] and sha(raw)==fc['source_txt_sha256']
  lines=raw.decode('ascii').splitlines();b,i,nf,nstate=map(int,lines[0].split());assert (b,i,nf)==(B,I,26-B) and len(lines)==nstate+1
  lp={tuple(map(int,k.split(','))) if k else ():v for k,v in fc['certified_signatures'].items()}
  for key,nums in lp.items():
   assert list(key)==sorted(set(key)) and all(0<m<(1<<B) for m in key)
   gap=verify_lp(B,key,nums);minlp=gap if minlp is None else min(minlp,gap)
  counts=Counter();surv=set();used=set();indices=[]
  for idx,line in enumerate(lines[1:]):
   ids=list(map(int,line.split(',')));assert len(ids)==nf and len(set(ids))==nf and all(0<=t<len(triples) for t in ids)
   faces=[list(triples[k]) for k in ids];key=signature(B,faces)
   if key is None:counts['direct']+=1
   elif key in lp:counts['farkas']+=1;used.add(key)
   else:
    counts['residual']+=1;surv.add(key);indices.append(idx);regenerated.append({'B':B,'I':I,'idx':idx,'faces':faces,'signature':list(key)})
  assert used==set(lp) and {','.join(map(str,k)) for k in surv}==set(fc['surviving_signatures'])
  assert dict(counts)==fc['counts'] and indices==fc['residual_indices'] and sum(counts.values())==nstate
  v=triangulation_count(B,I);num=2*factorial(2*B-3)*factorial(4*I+2*B-5);den=factorial(B-1)*factorial(B-3)*factorial(3*I+2*B-3)
  assert num%den==0 and v==num//den
  report[name]={'orbits':nstate,**counts,'labeled_count_by_root_peeling':v,'source_txt_sha256':sha(raw)}
  print(name,report[name],flush=True)
 assert regenerated==given
 candidates=[i for i,r in enumerate(given) if r['B']==10 and r['I']==4 and r['faces']==[list(f) for f in ALL_FACES]]
 assert candidates==[685]
 total=sum(r['orbits'] for r in report.values());remaining=len(given)
 out={'verified':True,'families':report,'total_orbits':total,'metric_excluded':total-remaining,'residuals':remaining,'candidate_residual_index':685,'candidate_enumeration_index':given[685]['idx'],'constants':{'A':A,'C':C,'D':D,'RU':str(RU)},'sin_margins':{k:str(v-(2*RU if k=='four_rod' else RU)) for k,v in sins.items()},'min_Farkas_margin':str(F(minlp,D*S)),'residual_sha256':sha((ENUM/'metric_residuals.json').read_bytes()),'seconds':time.time()-st}
 (BASE/'partition_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':verify()
