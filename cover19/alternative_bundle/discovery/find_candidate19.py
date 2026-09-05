from collections import deque
from pathlib import Path
import geometry19 as S,json
R=S.R

def canon(B,faces):
 inc={}
 for k,f in enumerate(faces):
  for i in range(3):inc.setdefault(tuple(sorted((f[i],f[(i+1)%3]))),[]).append(k)
 best=None;bestmp=None
 for dr in (1,-1):
  for start in range(B):
   mp={(start+dr*i)%B:i for i in range(B)};u,v=start,(start+dr)%B
   root=inc[tuple(sorted((u,v)))][0];todo=deque([(root,u,v)]);seen={root};out=[]
   while todo:
    j,a,b=todo.popleft();w=next(x for x in faces[j] if x not in (a,b))
    if w not in mp:mp[w]=len(mp)
    out.append(sum(1<<mp[x] for x in faces[j]))
    for edge in ((b,w),(w,a)):
     for jj in inc[tuple(sorted(edge))]:
      if jj not in seen:seen.add(jj);todo.append((jj,*edge))
   out=tuple(sorted(out))
   if best is None or out<best:best=out;bestmp=mp
 return best,bestmp

def profile(fs):
 adj=[set() for _ in range(19)]
 for f in fs:
  for v in f:adj[v].update(set(f)-{v})
 return sorted(map(len,adj[:12])),sorted(map(len,adj[12:]))
code,mp=canon(12,S.FACES);prof=profile(S.FACES);matches=[]
rows=json.load(open(R/'metric19_residuals.json'))
for k,row in enumerate(rows):
 if row['B']!=12 or profile(row['faces'])!=prof:continue
 co,m=canon(12,row['faces'])
 if co==code:
  inv={v:u for u,v in m.items()};p=[inv[mp[i]] for i in range(19)]
  matches.append({'residual_index':k,'map':p,'row':row})
assert len(matches)==1
(R/'candidate19_map.json').write_text(json.dumps(matches[0],indent=2));print(matches[0])
