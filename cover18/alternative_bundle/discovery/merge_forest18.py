from pathlib import Path
from collections import Counter
from fractions import Fraction
import json,gzip,time,hashlib
R=Path(__file__).resolve().parent;rows=json.load(open(R/'metric18_residuals.json'));cand=json.load(open(R/'candidate18_map.json'))['residual_index'];counts=Counter();maxdepth=0;st=time.time();mar=None
head={'version':1,'residual_sha256':hashlib.sha256((R/'metric18_residuals.json').read_bytes()).hexdigest(),'radius_upper':str(Fraction('0.29016771764058')),'cases':len(rows)-1}
with gzip.open(R/'noncandidate18_trees.jsonl.gz','wt',compresslevel=6,encoding='utf8') as f:
 f.write(json.dumps(head,separators=(',',':'))+'\n')
 for k,row in enumerate(rows):
  if k==cand:continue
  with gzip.open(R/'forest18'/f'{k:05}.json.gz','rt',encoding='utf8') as g:c=json.load(g)
  assert c['residual_index']==k and all(c[key]==row[key] for key in row)
  counts.update(n['kind'] for n in c['nodes']);maxdepth=max(maxdepth,c['summary']['maxdepth']);m=Fraction(c['summary']['min_margin']);mar=m if mar is None else min(mar,m)
  # Remove discovery-only optimizer summaries. The replay recomputes them.
  c.pop('summary',None);f.write(json.dumps(c,separators=(',',':'))+'\n')
  if k%5000==0:print('MERGE',k,'seconds',time.time()-st,flush=True)
print('MERGED',dict(counts),'maxdepth',maxdepth,'minmargin',float(mar),'seconds',time.time()-st,flush=True)
