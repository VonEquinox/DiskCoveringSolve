"""Untrusted discovery of exactly checkable angle/Farkas obstructions, n=15."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,math
import numpy as np
from numba import njit
from scipy.optimize import linprog
R=Path(__file__).resolve().parent
A=10306;C=21954;D=100000;S=10**12
RU=F('0.31814293085928')
@njit(cache=True)
def signatures(B,arr):
    n,nf=arr.shape;out=np.zeros((n,B*B),np.int32); counts=np.zeros(n,np.int32)
    for it in range(n):
        adj=np.array([1<<i for i in range(B)],np.int64)
        for z in arr[it]:
            mask=int(z)&((1<<B)-1)
            for u in range(B):
                if (mask>>u)&1:adj[u]|=mask
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
                if not present:out[it,kcount]=mask;kcount+=1
            if bad:break
        if bad:counts[it]=-1
        else:
            counts[it]=kcount;out[it,:kcount]=np.sort(out[it,:kcount])
    return out,counts

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def make():
 st=time.time();cert={};res=[]
 for B in range(10,15):
    p=R/f'B{B}_I{15-B}.masks';head=np.fromfile(p,dtype=np.uint32,count=5);assert int(head[0])==0x4449534b and int(head[1])==B and int(head[2])==15-B
    nf=int(head[3]);n=int(head[4]);arr=np.memmap(p,dtype=np.uint16,mode='r',offset=20,shape=(n,nf));assert nf==28-B
    direct=0;cache={};surv={};inds=[];nfark=0
    for start in range(0,n,50000):
     out,count=signatures(B,arr[start:start+50000]);direct+=int(np.sum(count<0))
     for off in np.flatnonzero(count>=0):
        idx=start+int(off);key=tuple(map(int,out[off,:count[off]]))
        if key not in cache:
            masks=[1<<i for i in range(B)]+list(key);mat=np.array([[bool(m>>i&1) for m in masks] for i in range(B)],float);cost=np.r_[np.full(B,A/D),np.full(len(key),C/D)]
            lp=linprog(cost,A_ub=-mat,b_ub=-np.ones(B),bounds=(0,None),method='highs');assert lp.success
            nums=[max(0,int(math.ceil(float(v)*S-1e-7))) for v in lp.x]
            for i in range(B):
                gap=S-sum(v for m,v in zip(masks,nums) if m>>i&1)
                if gap>0:nums[i]+=gap
            cc=A*sum(nums[:B])+C*sum(nums[B:]);cache[key]=nums if cc<D*S else None
        if cache[key] is not None:nfark+=1
        else:
            surv[key]=1;inds.append(idx);faces=[[i for i in range(15) if (int(z)>>i)&1] for z in arr[idx]]
            assert all(len(f)==3 for f in faces)
            res.append({'B':B,'I':15-B,'idx':idx,'faces':faces,'signature':list(key)})
     if start%1000000==0:print('progress',B,start,'direct',direct,'res',len(inds),'sec',round(time.time()-st,2),flush=True)
    family=f'B{B}_I{15-B}'
    cert[family]={'source_masks_sha256':sha(p),'counts':{'total':n,'direct':direct,'farkas':nfark,'residual':len(inds)},'certified_signatures':{','.join(map(str,k)):v for k,v in cache.items() if v is not None},'surviving_signatures':[','.join(map(str,k)) for k in surv],'residual_indices':inds}
    print(family,'total',n,'direct',direct,'Farkas',nfark,'residual',len(inds),'LPs',len(cache),'seconds',time.time()-st,flush=True)
    (R/'metric15_certificate.json').write_text(json.dumps(cert,separators=(',',':')));(R/'metric15_residuals.json').write_text(json.dumps(res,separators=(',',':')))
 print('COMPLETE',len(res),'seconds',time.time()-st,flush=True)
if __name__=='__main__':make()
