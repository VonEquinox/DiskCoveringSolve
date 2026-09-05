from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay
from scipy.optimize import minimize
import json, pathlib,sys,math,time
B=12;N=16

def data(points):
    fa=sorted(tuple(sorted(map(int,t))) for t in Delaunay(points).simplices)
    assert len(fa)==2*N-B-2
    e=[(i,B+i) for i in range(B)]+[(i,B+(i+1)%B) for i in range(B)]+[(B+N+j,B+i) for j,f in enumerate(fa) for i in f]
    return fa,np.array(e)

def solve(P,fa,E):
    # q_i is between c_i and c_(i+1); rotate q0 to (1,0)
    ag=np.unwrap(np.arctan2(P[:B,1],P[:B,0])); a0=(ag[0]+ag[1])/2
    rot=np.array([[np.cos(a0),-np.sin(a0)],[np.sin(a0),np.cos(a0)]])
    P=P@rot
    theta=np.r_[np.arange(1,B)*2*np.pi/B]
    W=[]
    for f in fa:
        a,b,c=P[list(f)];w=np.linalg.solve(2*np.array([b-a,c-a]),np.array([b@b-a@a,c@c-a@a]));W.append(w)
    x=np.r_[theta,P.ravel(),np.array(W).ravel(),.12]
    M=len(E);nn=B+N+len(fa); nv=len(x)
    def calc(x,jac=False):
        q=np.c_[np.cos(np.r_[0,x[:B-1]]),np.sin(np.r_[0,x[:B-1]])]
        z=np.r_[q,x[B-1:-1].reshape(-1,2)]
        dif=z[E[:,0]]-z[E[:,1]];g=x[-1]-np.sum(dif*dif,axis=1)
        if not jac:return g
        J=np.zeros((M,nv));J[:,-1]=1
        for k,(u,v) in enumerate(E):
            for w,sg in [(u,-1),(v,1)]:
                if w==0:continue
                if w<B:J[k,w-1]+=sg*2*dif[k]@np.array([-q[w,1],q[w,0]])
                else:J[k,B-1+2*(w-B):B-1+2*(w-B)+2]+=sg*2*dif[k]
        return J
    obj=np.zeros(nv);obj[-1]=1
    con={'type':'ineq','fun':calc,'jac':lambda x:calc(x,True)}
    res=minimize(lambda x:x[-1],x,jac=lambda x:obj,constraints=con,method='SLSQP',options={'ftol':1e-14,'maxiter':600,'disp':False})
    q=np.c_[np.cos(np.r_[0,res.x[:B-1]]),np.sin(np.r_[0,res.x[:B-1]])]
    z=np.r_[q,res.x[B-1:-1].reshape(-1,2)];lens=((z[E[:,0]]-z[E[:,1]])**2).sum(axis=1)
    return {'success':bool(res.success),'message':str(res.message),'t':float(res.fun),'r':float(np.sqrt(res.fun)),'theta':np.r_[0,res.x[:B-1]].tolist(),'points':z[B:B+N].tolist(),'witnesses':z[B+N:].tolist(),'faces':fa,'edges':E.tolist(),'rod_slacks':(res.fun-lens).tolist(),'iterations':int(res.nit)}

if __name__=='__main__':
    start=time.time();rng=np.random.default_rng(15)
    best=1.;path=pathlib.Path(__file__).parent
    for k in range(int(sys.argv[1]) if len(sys.argv)>1 else 10):
        th=np.arange(B)*2*np.pi/B
        P=np.r_[.82*np.c_[np.cos(th),np.sin(th)],.35*np.c_[np.cos(np.arange(4)*2*np.pi/4+.1*(k%5)),np.sin(np.arange(4)*2*np.pi/4+.1*(k%5))]]
        if k>=5:P+=rng.normal(0,.02,P.shape)
        fa,E=data(P)
        try:d=solve(P,fa,E)
        except Exception as ex:print(k,ex,flush=True);continue
        print(k,d['r'],d['iterations'],d['success'],'sec',time.time()-start,flush=True)
        (path/f'trial11_{k:03}.json').write_text(json.dumps(d,indent=2))
        if d['t']<best and min(d['rod_slacks'])>-1e-9:
            best=d['t'];(path/'candidate16_B12_numeric.json').write_text(json.dumps(d,indent=2))
