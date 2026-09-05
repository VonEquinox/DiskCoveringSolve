"""Stream every noncandidate certificate; independently reconstruct its domain.

The worker arithmetic is exact. Multiprocessing is only a deterministic map;
the parent checks the complete ordered case set before accepting any result.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from collections import Counter,deque
import gzip,json,time,hashlib,os
from cases20 import Store,ROOT,ENUM
import system20 as S
from verify_tree20 import verify as tree_verify
if not __debug__:raise RuntimeError('Run without -O')
STORE=None

def init():
 global STORE
 STORE=Store()

def batch(records):
 stats={'cases':0,'nodes':0,'counts':Counter(),'maxdepth':0,'minimum_force_margin':None,'scales':Counter(),'by_B':{}}
 for rec in records:
  assert type(rec)is dict and set(rec)=={'B','idx','nodes'}
  B,j=rec['B'],rec['idx'];assert type(B)is int and type(j)is int
  out=tree_verify(B,STORE.faces(B,j),rec['nodes'],False)
  assert out['counts']['LOCAL']==0
  stats['cases']+=1;stats['nodes']+=out['nodes'];stats['counts'].update(out['counts']);stats['scales'].update(out['scales'])
  stats['maxdepth']=max(stats['maxdepth'],out['maxdepth']);m=out['minimum_force_margin']
  if m is not None:stats['minimum_force_margin']=m if stats['minimum_force_margin'] is None else min(m,stats['minimum_force_margin'])
  q=stats['by_B'].setdefault(B,{'cases':0,'nodes':0,'counts':Counter()});q['cases']+=1;q['nodes']+=out['nodes'];q['counts'].update(out['counts'])
 return stats

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def aggregate(s,t):
 for k in ['cases','nodes']:s[k]+=t[k]
 for k in ['counts','scales']:s[k].update(t[k])
 s['maxdepth']=max(s['maxdepth'],t['maxdepth']);m=t['minimum_force_margin']
 if m is not None:s['minimum_force_margin']=m if s['minimum_force_margin'] is None else min(m,s['minimum_force_margin'])
 for B,q in t['by_B'].items():
  a=s['by_B'].setdefault(B,{'cases':0,'nodes':0,'counts':Counter()})
  a['cases']+=q['cases'];a['nodes']+=q['nodes'];a['counts'].update(q['counts'])

def plain(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [plain(v) for v in x]
 return x

def consume_key(expected,rec):
 key=next(expected,None)
 assert key is not None and (rec['B'],rec['idx'])==key,'missing, duplicated, or reordered case'
 return key

def finish_keys(expected):
 assert next(expected,None) is None,'missing final cases'

def verify(workers=4):
 start=time.time();store=Store();mp=json.load(open(ENUM/'candidate20_map.json'));candidate=(mp['B'],mp['idx'])
 cp=ROOT/'core/candidate20_tree.json.gz'
 with gzip.open(cp,'rt') as f:ct=json.load(f)
 assert set(ct)=={'version','faces','nodes'} and ct['version']==1 and ct['faces']==[list(f) for f in S.FACES]
 cand=tree_verify(S.B,S.FACES,ct['nodes'],True)
 assert cand['counts']['LOCAL']>0 and cand['counts']['EMPTY']==0
 cand['certificate_sha256']=sha(cp)
 print('CANDIDATE DOMAIN EXACTLY VERIFIED',cand['nodes'],cand['counts'],flush=True)
 expected=(key for key in store.keys() if key!=candidate)
 path=ENUM/'noncandidate20_trees.jsonl.gz'
 stats={'cases':0,'nodes':0,'counts':Counter(),'maxdepth':0,'minimum_force_margin':None,'scales':Counter(),'by_B':{}}
 sent=0;chunks=0;pending=deque()
 with ProcessPoolExecutor(max_workers=workers,initializer=init) as ex,gzip.open(path,'rt',encoding='utf8') as f:
  records=[]
  for line in f:
   assert line.strip(),'empty record'
   rec=json.loads(line);assert type(rec)is dict and set(rec)=={'B','idx','nodes'}
   consume_key(expected,rec)
   records.append(rec);sent+=1
   if len(records)==400:
    pending.append(ex.submit(batch,records));records=[]
    if len(pending)>=2*workers:
     aggregate(stats,pending.popleft().result());chunks+=1
     if chunks%100==0:print('NONCANDIDATE VERIFIED',stats['cases'],'/',sum(store.counts.values())-1,'nodes',stats['nodes'],'sec',round(time.time()-start,2),flush=True)
  finish_keys(expected)
  if records:pending.append(ex.submit(batch,records))
  while pending:aggregate(stats,pending.popleft().result())
 assert sent==stats['cases']==sum(store.counts.values())-1
 assert stats['counts']['LOCAL']==0 and stats['counts']['EMPTY']==0
 assert stats['counts']['DUAL']==stats['counts']['SPLIT']+stats['cases']
 assert stats['minimum_force_margin']>0
 stats['certificate_sha256']=sha(path)
 out={'verified':True,'N':20,'candidate_key':list(candidate),'candidate_tree':cand,'noncandidate_forest':stats,
      'surviving_orbits':sum(store.counts.values()),'unresolved_leaves':0,
      'root_sha256':sha(ROOT/'core/root20.json'),'local_certificate_sha256':sha(ROOT/'core/anchor_isolation_certificate.json'),
      'seconds':time.time()-start}
 (ROOT/'reports/forest20_verified.json').write_text(json.dumps(plain(out),indent=2))
 print('ALL CONTINUOUS DOMAINS EXACTLY VERIFIED',stats['cases'],stats['nodes'],dict(stats['counts']),'candidate',cand['nodes'],'seconds',round(time.time()-start,2),flush=True)
 return out
if __name__=='__main__':verify(int(os.environ.get('COVER20_WORKERS','4')))
