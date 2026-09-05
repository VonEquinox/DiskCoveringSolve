from pathlib import Path
import json
R=Path(__file__).resolve().parent
rows=json.load(open(R/'metric15_residuals.json'))
for b in range(10,15):
 (R/f'metric15_residuals_B{b}.json').write_text(json.dumps([r for r in rows if r['B']==b],separators=(',',':')))
