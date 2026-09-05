"""Create a separate optional discovery workspace; never modify proof inputs."""
from pathlib import Path
import os,shutil
P=Path(__file__).resolve().parents[1];W=P/'_discovery_work';W.mkdir(exist_ok=True)
for source in [P/'core',P/'discovery']:
 for p in source.iterdir():
  if p.is_file() and p.suffix in ['.py','.json'] and p.name!='prepare_workspace.py':shutil.copy2(p,W/p.name)
for name in ['metric15_residuals.json','metric15_certificate.json','candidate15_map.json']:shutil.copy2(P/'enumeration'/name,W/name)
for B in range(10,15):
 src=P/'_runtime'/f'B{B}_I{15-B}.masks';dst=W/src.name
 if not src.is_file():raise FileNotFoundError('Run verify_all.py first: '+str(src))
 if dst.exists():dst.unlink()
 try:os.link(src,dst)
 except OSError:shutil.copy2(src,dst)
print(W)
