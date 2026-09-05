from pathlib import Path
import numpy as np,json
from scipy.linalg import null_space
p=json.load(open(str(Path(__file__).resolve().parent)+'/candidate14_numeric.json'));P=np.array(p['points']);fa=p['faces'];active=[2,3,5,7,9,10,12,13]
q=np.c_[np.cos(p['theta']),np.sin(p['theta'])];W=[]
for i in active:
 a,b,c=P[fa[i]];W.append(np.linalg.solve(2*np.array([b-a,c-a]),np.array([b@b-a@a,c@c-a@a])))
z=np.r_[q,P,W];ed=[(i,10+i) for i in range(10)]+[(i,10+(i+1)%10) for i in range(10)]+[(24+j,10+i) for j,k in enumerate(active) for i in fa[k]]
dim=2*(len(z)-1)+1;J=[]
for u,v in ed:
 row=np.zeros(dim)
 for a,b in [(u,v),(v,u)]:
  if a:row[2*(a-1):2*a]=2*(z[a]-z[b])
 row[-1]=-1;J.append(row)
for i in range(1,10):
 row=np.zeros(dim);row[2*(i-1):2*i]=2*z[i];J.append(row)
J=np.array(J);e=np.zeros(dim);e[-1]=-1
w=np.linalg.lstsq(J.T,e,rcond=1e-12)[0]
L=np.zeros((len(z),len(z)))
for (u,v),a in zip(ed,w):L[u,u]+=a;L[v,v]+=a;L[u,v]-=a;L[v,u]-=a
for i in range(1,10):L[i,i]+=w[len(ed)+i-1]
M=np.kron(L[1:,1:],np.eye(2)); A=J[:,:-1]
K=null_space(A,rcond=1e-7)
print('active',active,'rank J A',np.linalg.matrix_rank(J),np.linalg.matrix_rank(A,tol=1e-7),'res',max(abs(J.T@w-e)))
print('last singular',np.linalg.svd(A)[1][-5:]);print('tangent eig',np.linalg.eigvalsh(K.T@M@K))
print('min w',min(w[:len(ed)]),'mu',w[len(ed):]);print('Mnorm',np.linalg.eigvalsh(M)[[0,-1]])
pa=np.zeros_like(M);pa[:18,:18]=np.eye(18)
wa=np.r_[w[:len(ed)]**2,np.full(9,.02**2)]
for c in [10,100,150,500,1000,10000]:
 print('penalty',c, np.linalg.eigvalsh(M+c*A.T@np.diag(wa)@A-.001*pa)[:4])
out={'B':10,'N':14,'active_faces':[fa[i] for i in active],'all_faces':fa,'edges':ed,'x':np.r_[z[1:].ravel(),p['t'],w].tolist(),'geometry_count':dim-1,'constraints':len(J)}
open(str(Path(__file__).resolve().parent)+'/active14_numeric.json','w').write(json.dumps(out,indent=2))
