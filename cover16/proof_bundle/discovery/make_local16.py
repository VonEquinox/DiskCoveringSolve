from pathlib import Path
import numpy as np,json,hashlib
import anchor_isolation16 as A
H,err,cert=A.prepare();z=np.array([[float(x/A.DEN) for x in row] for row in H]);print('min eigen',np.linalg.eigvalsh(z)[0]);L=np.linalg.cholesky(z);R=np.linalg.inv(L.T);RN=[[int(round(v*A.QR)) if j>=i else 0 for j,v in enumerate(row)] for i,row in enumerate(R)]
c={'QR':A.QR,'Rnum':RN,'root_sha256':A.h(A.ROOT/'root16.json'),'root_data_sha256':A.h(A.ROOT/'system16.py')}
(A.ROOT/'anchor_isolation_certificate.json').write_text(json.dumps(c,separators=(',',':')))
A.main()
