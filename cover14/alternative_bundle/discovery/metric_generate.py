"""Numerical LP discovery with independent exact checks of every accepted LP certificate."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
import gzip,json,hashlib,time,numpy as np
from scipy.optimize import linprog
BASE=Path(__file__).resolve().parent;ENUM=BASE/'enumeration'
A=10764;C=23092;D=100000;S=10**12;N=14

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
   if fl:mask=((1<<j)-1)^((1<<i)-1)
   elif fr:mask=((1<<B)-1)^(((1<<j)-1)^((1<<i)-1))
   else:continue
   masks.add(mask)
 return tuple(sorted(masks))

def verify_lp(B,key,nums):
 masks=[1<<i for i in range(B)]+list(key)
 assert len(nums)==len(masks) and all(type(v)is int and v>=0 for v in nums)
 assert all(sum(v for v,m in zip(nums,masks) if m>>i&1)>=S for i in range(B))
 cost=A*sum(nums[:B])+C*sum(nums[B:]);assert cost<D*S
 return D*S-cost

def lp(B,key):
 masks=[1<<i for i in range(B)]+list(key)
 mat=np.array([[(mask>>i)&1 for mask in masks] for i in range(B)],float)
 rhs=np.r_[np.full(B,A),np.full(len(key),C)]/D
 out=linprog(rhs,A_ub=-mat,b_ub=-np.ones(B),bounds=(0,None),method='highs')
 if out.success and out.fun<1-1e-9:
  nums=[int(np.ceil(v*S+1e-2)) for v in out.x];verify_lp(B,key,nums);return nums
 return None

def main():
 st=time.time();total={};res=[];trips=list(combinations(range(14),3))
 for B in range(10,14):
  fam=f'B{B}_I{14-B}';raw=(ENUM/(fam+'.jsonl.txt')).read_bytes();packed=gzip.compress(raw,mtime=0)
  (ENUM/(fam+'.txt.gz')).write_bytes(packed)
  lines=raw.decode().splitlines();b,i,nf,nstates=map(int,lines[0].split());assert (b,i,nf)==(B,14-B,26-B)
  cache={};cnt=Counter();inds=[];allsig=Counter()
  for idx,line in enumerate(lines[1:]):
   faces=[list(trips[int(s)]) for s in line.split(',')];key=signature(B,faces)
   if key is None:cnt['direct']+=1;continue
   if key not in cache:cache[key]=lp(B,key)
   if cache[key] is not None:cnt['farkas']+=1
   else:
    cnt['residual']+=1;inds.append(idx);res.append({'B':B,'I':14-B,'idx':idx,'faces':faces,'signature':list(key)});allsig[key]+=1
  total[fam]={'counts':dict(cnt),'certified_signatures':{','.join(map(str,key)):val for key,val in cache.items() if val is not None},'surviving_signatures':[','.join(map(str,key)) for key,val in cache.items() if val is None],'residual_indices':inds,'source_txt_sha256':hashlib.sha256(raw).hexdigest(),'source_gzip_sha256':hashlib.sha256(packed).hexdigest()}
  print(fam,cnt,'signatures',len(cache),'surviving',len(allsig),'elapsed',round(time.time()-st,1),flush=True)
  (ENUM/'metric_partition_certificate.json').write_text(json.dumps(total,separators=(',',':')))
  (ENUM/'metric_residuals.json').write_text(json.dumps(res,separators=(',',':')))
 print('DONE',len(res),flush=True)
if __name__=='__main__':main()
