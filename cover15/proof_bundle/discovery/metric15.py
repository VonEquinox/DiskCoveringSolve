from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,math,sys
import numpy as np
from numba import njit
from scipy.optimize import linprog
R=Path(__file__).resolve().parent
A=10306;C=21954;D=100000;S=10**12
RU=F('0.318142930859264')
@njit
def signatures(B,arr):
    n,nf=arr.shape;out=np.zeros((n,2*B),np.int64); counts=np.zeros(n,np.int64)
    for it in range(n):
        adj=np.array([1<<i for i in range(B)],np.int64)
        for z in arr[it]:
            a=z//225;b=(z//15)%15;c=z%15
            t=np.array([a,b,c])
            for u in t:
                if u<B:
                    for v in t:
                        if v<B:adj[u]|=1<<v
        kcount=0;bad=False
        for i in range(B):
            reach=adj[i]|adj[(i+1)%B]
            for j in range(i+1,B):
                if not(reach&((1<<j)|(1<<((j+1)%B)))):continue
                k=j-i
                if min(k,B-k)<=2:continue
                fl=k*A<D-C;fr=(B-k)*A<D-C
                if fl and fr:bad=True;break
                if fl:mask=((1<<j)-1)^((1<<i)-1)
                elif fr:mask=((1<<B)-1)^(((1<<j)-1)^((1<<i)-1))
                else:continue
                present=False
                for m in range(kcount):
                    if out[it,m]==mask:present=True;break
                if not present:
                    assert kcount<2*B
                    out[it,kcount]=mask;kcount+=1
            if bad:break
        if bad:counts[it]=-1
        else:
            counts[it]=kcount;out[it,:kcount]=np.sort(out[it,:kcount])
    return out,counts

def make():
 st=time.time();cert={};res=[]
 for B in ([int(sys.argv[1])] if len(sys.argv)>1 else range(10,15)):
    p=R/f'B{B}_I{15-B}.bin';head=np.fromfile(p,dtype=np.uint32,count=5);nf=int(head[3]);n=int(head[4]);arr=np.fromfile(p,dtype=np.uint16,offset=20).reshape(n,nf)
    out,count=signatures(B,arr);direct=int(np.sum(count<0));cache={};surv={};inds=[];nfark=0
    for idx in np.flatnonzero(count>=0):
        key=tuple(map(int,out[idx,:count[idx]]))
        if key not in cache:
            masks=[1<<i for i in range(B)]+list(key);mat=np.array([[bool(m>>i&1) for m in masks] for i in range(B)],float);cost=np.r_[np.full(B,A/D),np.full(len(key),C/D)]
            lp=linprog(cost,A_ub=-mat,b_ub=-np.ones(B),bounds=(0,None),method='highs')
            assert lp.success
            nums=np.ceil(np.maximum(lp.x,0)*S-1e-7).astype(np.int64)
            # Ensure exact coverage after rounding (integer correction).
            for i in range(B):
                gap=S-sum(int(v) for m,v in zip(masks,nums) if m>>i&1)
                if gap>0:nums[i]+=gap
            cc=A*sum(map(int,nums[:B]))+C*sum(map(int,nums[B:]))
            cache[key]=nums.tolist() if cc<D*S else None
        if cache[key] is not None:nfark+=1
        else:
            surv[key]=1;inds.append(int(idx));faces=[[int(z)//225,(int(z)//15)%15,int(z)%15] for z in arr[idx]]
            res.append({'B':B,'I':15-B,'idx':int(idx),'faces':faces,'signature':list(key)})
    family=f'B{B}_I{15-B}';txt=R/(family+'.txt')
    cert[family]={'source_txt_sha256':hashlib.sha256(txt.read_bytes()).hexdigest(),'counts':{'direct':direct,'farkas':nfark,'residual':len(inds)},'certified_signatures':{','.join(map(str,k)):v for k,v in cache.items() if v is not None},'surviving_signatures':[','.join(map(str,k)) for k in surv],'residual_indices':inds}
    print(family,'total',n,'direct',direct,'Farkas',nfark,'residual',len(inds),'LPs',len(cache),'seconds',time.time()-st,flush=True)
 (R/(f'metric15_certificate_B{sys.argv[1]}.json' if len(sys.argv)>1 else 'metric15_certificate.json')).write_text(json.dumps(cert,separators=(',',':')));(R/(f'metric15_residuals_B{sys.argv[1]}.json' if len(sys.argv)>1 else 'metric15_residuals.json')).write_text(json.dumps(res,separators=(',',':')))
if __name__=='__main__':make()
