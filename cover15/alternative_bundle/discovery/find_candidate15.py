from pathlib import Path
import system15 as S,json,itertools
rows=json.load(open(S.R/'metric15_residuals.json'));dic={tuple(sorted(tuple(f) for f in row['faces'])):k for k,row in enumerate(rows) if row['B']==11};matches=[]
for rev in [False,True]:
 for sh in range(11):
  for ip in itertools.permutations(range(11,15)):
   mp=[((sh-i if rev else sh+i)%11) for i in range(11)]+list(ip);faces=tuple(sorted(tuple(sorted(mp[v] for v in f)) for f in S.FACES))
   if faces in dic:matches.append({'residual_index':dic[faces],'map':mp,'reflection':rev,'shift':sh})
assert matches and len(set(r['residual_index'] for r in matches))==1
out=matches[0];out['row']=rows[out['residual_index']];print(out)
(S.R/'candidate15_map.json').write_text(json.dumps(out,indent=2))
