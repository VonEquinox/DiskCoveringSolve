import mpmath as mp,json,time
from pathlib import Path
import system16 as S
mp.mp.dps=110
R=Path(__file__).resolve().parent
d=json.load(open(R/'active16_numeric.json'));z=[mp.mpf(str(x)) for x in d['x']]
t=time.time()
for k in range(8):
    f,j=S.evaluate(z,mp.mpf(0));err=max(map(abs,f));print(k,mp.nstr(err,12),flush=True)
    if err<mp.mpf('1e-100'):break
    step=mp.lu_solve(mp.matrix(j),mp.matrix(f));z=[x-dx for x,dx in zip(z,step)]
f,j=S.evaluate(z,mp.mpf(0));Y=mp.inverse(mp.matrix(j));qx=10**90;qy=10**80
r=mp.sqrt(z[S.TID]);print('r',mp.nstr(r,105),'R',mp.nstr(1/r,105),flush=True)
out={'version':1,'N':S.D,'Qx':str(qx),'Qy':str(qy),'rho':'1/'+str(10**60),'xnum':[str(int(mp.nint(x*qx))) for x in z],'Ynum':[[str(int(mp.nint(Y[i,j]*qy))) for j in range(S.D)] for i in range(S.D)],'r_decimal':mp.nstr(r,100),'R_decimal':mp.nstr(1/r,100),'active_faces':S.ACTIVE,'all_faces':S.FACES,'edges':S.EDGES}
(R/'root16.json').write_text(json.dumps(out,separators=(',',':')))
P=[mp.matrix([mp.mpf(1),mp.mpf(0)])]+[mp.matrix(z[2*i:2*i+2]) for i in range(S.G//2)]
for f in S.FACES:
    a,b,c=[P[S.B+i] for i in f];A=mp.matrix([[2*(b-a)[i] for i in (0,1)],[2*(c-a)[i] for i in (0,1)]]);rhs=mp.matrix([(b.T*b)[0]-(a.T*a)[0],(c.T*c)[0]-(a.T*a)[0]])
    q=mp.lu_solve(A,rhs);rr=((q-a).T*(q-a))[0];print(f,mp.nstr(z[S.TID]-rr,60),flush=True)
print('elapsed',time.time()-t,flush=True)
