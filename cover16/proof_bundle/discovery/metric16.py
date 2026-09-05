"""Untrusted discovery of exactly checkable angle/Farkas obstructions, n=16."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,math
import numpy as np
from numba import njit
from scipy.optimize import linprog
R=Path(__file__).resolve().parent
A=9974;C=21143;C6=37566;D=100000;S=10**12
RU=F('0.30821980189863')
@njit(cache=True)
def signatures(B,arr):
    n,nf=arr.shape;out=np.zeros((n,B*B),np.int32); counts=np.zeros(n,np.int32)
    for it in range(n):
        adj=np.array([1<<i for i in range(16)],np.int64)
        for z in arr[it]:
            mask=int(z)
            for u in range(16):
                if (mask>>u)&1:adj[u]|=mask
        r2=adj.copy()
        for u in range(16):
            for v in range(16):
                if adj[u]>>v&1:r2[u]|=adj[v]
        kcount=0;bad=False
        for i in range(B):
            reach=adj[i]|adj[(i+1)%B];reach2=r2[i]|r2[(i+1)%B]
            for j in range(i+1,B):
                ends=(1<<j)|(1<<((j+1)%B))
                if reach&ends:cap=C;tag=0
                elif reach2&ends:cap=C6;tag=1<<B
                else:continue
                k=j-i
                if min(k,B-k)*A<=cap:continue
                fl=k*A<D-cap;fr=(B-k)*A<D-cap
                if fl and fr:bad=True;break
                if fl:mask=((1<<j)-1)^((1<<i)-1)
                elif fr:mask=((1<<B)-1)^(((1<<j)-1)^((1<<i)-1))
                else:continue
                mask|=tag
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
 for B in range(15,10,-1):
    while not (R/f'B{B}_I{16-B}.json').is_file():time.sleep(3)
    p=R/f'B{B}_I{16-B}.masks';head=np.fromfile(p,dtype=np.uint32,count=5);assert int(head[0])==0x4449534b and int(head[1])==B and int(head[2])==16-B
    nf=int(head[3]);n=int(head[4]);arr=np.memmap(p,dtype=np.uint16,mode='r',offset=20,shape=(n,nf));assert nf==30-B
    direct=0;cache={};surv={};inds=[];nfark=0
    for start in range(0,n,50000):
     out,count=signatures(B,arr[start:start+50000]);direct+=int(np.sum(count<0))
     for off in np.flatnonzero(count>=0):
        idx=start+int(off);key=tuple(map(int,out[off,:count[off]]))
        if key not in cache:
            from cycle_metric16 import cycle_weights
            ww=cycle_weights(B,np.array(key,dtype=np.int64))
            cache[key]=list(map(int,ww)) if len(ww) else None
        if cache[key] is not None:nfark+=1
        else:
            surv[key]=1;inds.append(idx);faces=[[i for i in range(16) if (int(z)>>i)&1] for z in arr[idx]]
            assert all(len(f)==3 for f in faces)
            res.append({'B':B,'I':16-B,'idx':idx,'faces':faces,'signature':list(key)})
     if start%1000000==0:print('progress',B,start,'direct',direct,'res',len(inds),'sec',round(time.time()-st,2),flush=True)
    family=f'B{B}_I{16-B}'
    cert[family]={'source_masks_sha256':sha(p),'counts':{'total':n,'direct':direct,'farkas':nfark,'residual':len(inds)},'certified_signatures':{','.join(map(str,k)):v for k,v in cache.items() if v is not None},'surviving_signatures':[','.join(map(str,k)) for k in surv],'residual_indices':inds}
    print(family,'total',n,'direct',direct,'Farkas',nfark,'residual',len(inds),'LPs',len(cache),'seconds',time.time()-st,flush=True)
    (R/'metric16_certificate.json').write_text(json.dumps(cert,separators=(',',':')));(R/'metric16_residuals.json').write_text(json.dumps(res,separators=(',',':')))
 print('COMPLETE',len(res),'seconds',time.time()-st,flush=True)
if __name__=='__main__':make()
