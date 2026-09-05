import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import anchor_isolation19 as A,numpy as np,json
H,den=A.prepare();Hf=np.array([[float(v/den) for v in row] for row in H]);eig=np.linalg.eigvalsh(Hf);print('min eigen',eig[:3])
T=np.linalg.inv(np.linalg.cholesky(Hf).T);nums=[[int(round(v*A.QR)) if j>=i else 0 for j,v in enumerate(row)] for i,row in enumerate(T)]
c={'QR':A.QR,'Rnum':nums,'geometry_sha256':A.sha(A.R/'geometry19.py')}
(A.R/'anchor_isolation_certificate.json').write_text(json.dumps(c,separators=(',',':')));A.verify()
