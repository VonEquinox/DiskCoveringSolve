#!/usr/bin/env python3
"""Recompute exact metric partition directly from the audited orbit representatives."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
import json,gzip,hashlib,time
if not __debug__:
 raise RuntimeError('Exact checks require Python without -O/PYTHONOPTIMIZE')
from exact_arcs import sincos_turn,SCALE
ROOT=Path(__file__).resolve().parent;ENUM=ROOT.parent/'enumeration'
S=10**12;A=1127;C=2440;D=10000
RU=F(346645456927389644,10**18)

def signature(B,faces):
 adj=[1<<i for i in range(B)]
 for fa in faces:
  bs=[v for v in fa if v<B]
  for u in bs:
   for v in bs:adj[u]|=1<<v
 masks=set()
 for i in range(B):
  reach=adj[i]|adj[(i+1)%B]
  for j in range(i+1,B):
   if not(reach&((1<<j)|(1<<((j+1)%B)))):continue
   k=j-i
   if min(k,B-k)<=2:continue
   fl=k*A<D-C;fr=(B-k)*A<D-C
   if fl and fr:return None
   if fl:mask=sum(1<<h for h in range(i,j))
   elif fr:mask=sum(1<<h for h in range(B) if not i<=h<j)
   else:continue
   masks.add(mask)
 return tuple(sorted(masks))

def verify_lp(B,key,nums):
 masks=[1<<i for i in range(B)]+list(key)
 assert len(nums)==len(masks) and all(type(v)is int and v>=0 for v in nums)
 assert all(sum(v for v,m in zip(nums,masks) if m>>i&1)>=S for i in range(B))
 cost=A*sum(nums[:B])+C*sum(nums[B:]);assert cost<D*S
 return D*S-cost

def check_constants():
 d=json.loads((ROOT/'cover13_krawczyk_cert.json').read_text());t=F(int(d['xnum'][62]),int(d['Qx']));rho=F(int(d['rho_num']),int(d['rho_den']))
 assert 0<rho and 0<t-rho<t+rho<RU**2<F(1,4)
 _,s1=sincos_turn(F(A,2*D));_,s2=sincos_turn(F(C,2*D));_,s8=sincos_turn(F(1,16))
 assert RU<F(s1[0],SCALE)
 assert 2*RU<F(s2[0],SCALE)
 assert RU<F(s8[0],SCALE)
 assert 0<2*A<C<D//2
 return {'radius_upper':str(RU),'cap_turns':str(F(A,D)),'four_rod_cap_turns':str(F(C,D)),
 'sine_cap_margin':str(F(s1[0],SCALE)-RU),'sine_four_rod_margin':str(F(s2[0],SCALE)-2*RU)}

def verify():
 st=time.time();constants=check_constants();cp=ENUM/'metric_partition_certificate.json';cert=json.loads(cp.read_text());triples=list(combinations(range(13),3))
 allres=[];summary={}
 for B in range(9,13):
  I=13-B;family=f'B{B}_I{I}';p=ENUM/(family+'.txt.gz');packed=p.read_bytes();raw=gzip.decompress(packed);c=cert[family]
  assert c['source_gzip_sha256']==hashlib.sha256(packed).hexdigest()
  assert c['source_txt_sha256']==hashlib.sha256(raw).hexdigest()
  lines=raw.decode('ascii').splitlines();b,i,nf,nstates=map(int,lines[0].split());assert (b,i,nf)==(B,I,24-B) and len(lines)-1==nstates
  cache={tuple(map(int,k.split(','))) if k else ():v for k,v in c['certified_signatures'].items()}
  surviving_expected={tuple(map(int,k.split(','))) if k else () for k in c['surviving_signatures']}
  for key,v in cache.items():verify_lp(B,key,v)
  counts=Counter();res=[];surviving=set()
  for idx,line in enumerate(lines[1:]):
   ids=list(map(int,line.split(',')));assert len(ids)==nf and len(set(ids))==nf and all(0<=a<len(triples) for a in ids)
   faces=[list(triples[a]) for a in ids];key=signature(B,faces)
   if key is None:counts['direct']+=1;continue
   if key in cache:counts['farkas']+=1
   else:
    assert key in surviving_expected
    counts['residual']+=1;surviving.add(key)
    res.append({'B':B,'I':I,'idx':idx,'faces':faces,'signature':list(key)})
  assert dict(counts)==c['counts'] and surviving==surviving_expected
  assert [r['idx'] for r in res]==c['residual_indices']
  allres.extend(res);summary[family]=dict(counts)
 assert allres==json.loads((ENUM/'metric_residuals_new.json').read_text())
 assert sum(sum(v.values()) for v in summary.values())==308198 and len(allres)==1761
 out={'verified':True,'constants':constants,'families':summary,'total':308198,'residual':1761,'certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'seconds':time.time()-st}
 (ROOT/'partition_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':verify()
