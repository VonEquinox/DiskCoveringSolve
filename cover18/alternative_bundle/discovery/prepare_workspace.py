"""Copy optional discovery tools into a fresh, separate flat working directory."""
from pathlib import Path
import argparse,shutil
ap=argparse.ArgumentParser();ap.add_argument('destination');args=ap.parse_args()
base=Path(__file__).resolve().parents[1];dst=Path(args.destination).expanduser().resolve()
if dst.exists():raise SystemExit('Destination must not exist; refusing to overwrite files')
if dst==base or base in dst.parents:raise SystemExit('Use a workspace outside the proof package')
dst.mkdir(parents=True)
for folder,patterns in [(base/'core',['*.py','root18.json','anchor_angles18.json','anchor_isolation_certificate.json']),
                        (base/'discovery',['*.py','*.json'])]:
    for pattern in patterns:
        for p in folder.glob(pattern):
            if p.name!='prepare_workspace.py':shutil.copyfile(p,dst/p.name)
shutil.copyfile(base/'enumeration'/'survivors18.json',dst/'metric18_residuals.json')
shutil.copyfile(base/'enumeration'/'candidate18_map.json',dst/'candidate18_map.json')
print('Created untrusted numerical-discovery workspace:',dst)
