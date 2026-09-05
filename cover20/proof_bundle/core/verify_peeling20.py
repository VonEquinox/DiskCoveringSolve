"""Compare two freshly run exhaustive recursions with all stored survivor data."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
from functools import lru_cache
import json,hashlib,subprocess,time
from cases20 import Store,ROOT,ENUM
import domain20 as D
if not __debug__:raise RuntimeError('Run without -O')
RUNTIME=ROOT/'_runtime';REPORTS=ROOT/'reports'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def labelled(B,I):
 assert B>=2 and I>=0
 if B==2:return int(I==0)
 z=sum(comb(I,j)*labelled(k,j)*labelled(B-k+1,I-j) for k in range(2,B) for j in range(I+1))
 if I:z+=I*(labelled(B+1,I-1)-sum(comb(I-1,j)*labelled(3,j)*labelled(B,I-1-j) for j in range(I)))
 return z

def verify():
 st=time.time();store=Store();summary=[];total=rooted=0
 for B in range(12,20):
  primary=json.load(open(RUNTIME/f'dfs_B{B}.json'));other=json.load(open(RUNTIME/f'bfs_B{B}.json'))
  for data in (primary,other):
   assert data['N']==20 and data['B']==B and data['I']==20-B and data['screen'] is True
   assert data['angle_scale']==D.D and data['angle_caps']==[D.A,D.C4,D.C6]
  assert other['verified'] is True
  rows=primary['survivors'];assert len(rows)==store.counts[B]==primary['distinct_survivors']==other['surviving_orbits']
  assert primary['complete_rooted']==other['rooted_survivors']
  dfs=RUNTIME/f'compare_dfs_B{B}.txt';bfs=RUNTIME/f'compare_bfs_B{B}.txt'
  with open(dfs,'w') as f:
   f.write(str(len(rows))+'\n')
   for j,ms in enumerate(rows):
    assert all(type(m)is int for m in ms) and tuple(ms)==store.masks(B,j)
    f.write(' '.join(map(str,ms))+'\n')
  with open(bfs,'w') as f:
   f.write(str(len(other['orbits']))+'\n')
   for row in other['orbits']:
    assert set(row)=={'masks','stabilizer','rooted_multiplicity'}
    assert all(type(m)is int for m in row['masks']) and type(row['stabilizer'])is int and type(row['rooted_multiplicity'])is int
    f.write(' '.join(map(str,[row['stabilizer'],row['rooted_multiplicity']]+row['masks']))+'\n')
  out=subprocess.run([str(RUNTIME/'compare_maps20'),'20',str(B),str(dfs),str(bfs)],capture_output=True,text=True,check=True)
  cmp=json.loads(out.stdout);assert cmp['verified'] is True and cmp['N']==20 and cmp['B']==B and cmp['orbits']==len(rows) and cmp['rooted']==primary['complete_rooted']
  lab=labelled(B,20-B);closed=2*factorial(2*B-3)*factorial(4*(20-B)+2*B-5)//(factorial(B-1)*factorial(B-3)*factorial(3*(20-B)+2*B-3));assert lab==closed
  summary.append({'B':B,'I':20-B,'orbits':len(rows),'rooted':cmp['rooted'],'dfs_calls':primary['calls'],'dfs_cuts':primary['prunes'],'bfs_calls':other['calls'],'bfs_cuts':other['metric_rejections'],'unpruned_labelled_count':lab,'survivors_sha256':sha(ENUM/f'survivors20_B{B}.bin.gz'),'dfs_sha256':sha(RUNTIME/f'dfs_B{B}.json'),'bfs_sha256':sha(RUNTIME/f'bfs_B{B}.json')})
  total+=len(rows);rooted+=cmp['rooted'];dfs.unlink();bfs.unlink()
  print('COMPLETE SURVIVOR SET VERIFIED',B,len(rows),cmp['rooted'],flush=True)
 out={'verified':True,'N':20,'families':summary,'surviving_orbits':total,'rooted_survivors':rooted,'unresolved_combinatorial_branches':0,'seconds':time.time()-st}
 (REPORTS/'peeling20_verified.json').write_text(json.dumps(out,indent=2));return out
if __name__=='__main__':verify()
