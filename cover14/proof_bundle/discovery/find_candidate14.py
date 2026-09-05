from pathlib import Path
import system14 as S,json,itertools
rows=json.load(open(S.R/'metric14_residuals.json'));dic={tuple(sorted(tuple(f) for f in row['faces'])):k for k,row in enumerate(rows) if row['B']==10};matches=[]
for rev in [False,True]:
 for sh in range(10):
  for ip in itertools.permutations(range(10,14)):
   mp=[((sh-i if rev else sh+i)%10) for i in range(10)]+list(ip);faces=tuple(sorted(tuple(sorted(mp[v] for v in f)) for f in S.FACES))
   if faces in dic:matches.append({'residual_index':dic[faces],'map':mp,'reflection':rev,'shift':sh})
assert matches and len(set(r['residual_index'] for r in matches))==1
out=matches[0];out['row']=rows[out['residual_index']];print(out)
(S.R/'candidate14_map.json').write_text(json.dumps(out,indent=2))
