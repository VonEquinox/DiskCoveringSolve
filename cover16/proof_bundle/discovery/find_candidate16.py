from pathlib import Path
import system16 as S,json,itertools
rows=json.load(open(S.R/'metric16_residuals.json'));dic={tuple(sorted(tuple(f) for f in row['faces'])):k for k,row in enumerate(rows) if row['B']==S.B};matches=[]
for rev in [False,True]:
 for sh in range(S.B):
  for ip in itertools.permutations(range(S.B,S.N)):
   mp=[((sh-i if rev else sh+i)%S.B) for i in range(S.B)]+list(ip);faces=tuple(sorted(tuple(sorted(mp[v] for v in f)) for f in S.FACES))
   if faces in dic:matches.append({'residual_index':dic[faces],'map':mp,'reflection':rev,'shift':sh})
assert matches and len(set(r['residual_index'] for r in matches))==1
out=matches[0];out['row']=rows[out['residual_index']];print(out)
(S.R/'candidate16_map.json').write_text(json.dumps(out,indent=2))
