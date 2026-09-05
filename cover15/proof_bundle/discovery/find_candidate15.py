from pathlib import Path
from itertools import permutations
import json
import system15 as S
R=Path(__file__).resolve().parent
rows=json.load(open(R/'metric15_residuals_B11.json'));idx={tuple(sorted(tuple(f) for f in row['faces'])):k for k,row in enumerate(rows)}
hits=[]
for rev in [False,True]:
 for sh in range(S.B):
  for ip in permutations(range(S.B,S.N)):
   mp=[((-i if rev else i)+sh)%S.B for i in range(S.B)]+list(ip)
   key=tuple(sorted(tuple(sorted(mp[v] for v in f)) for f in S.FACES))
   if key in idx:hits.append({'B':S.B,'I':S.N-S.B,'family_residual_index':idx[key],'enum_idx':rows[idx[key]]['idx'],'center_map':mp,'reflection':rev,'shift':sh})
assert hits and len(set(x['family_residual_index'] for x in hits))==1
(R/'candidate15_map_family.json').write_text(json.dumps(hits[0],indent=2))
print(hits)
