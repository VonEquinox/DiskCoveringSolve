"""Create a separate numerical discovery workspace without touching the proof.
Requires NumPy; subsequent proposal scripts also use SciPy/MPMath/Numba.
Usage: python bootstrap.py /path/to/new-empty-workspace
"""
from pathlib import Path
import sys,shutil,json,struct,gzip
import numpy as np
ROOT=Path(__file__).resolve().parents[1];SRC=ROOT/'discovery'
def main(path):
 out=Path(path).expanduser().resolve()
 if out==ROOT or ROOT in out.parents:raise ValueError('Choose a workspace outside the proof directory')
 out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Workspace must be empty')
 for p in SRC.iterdir():
  if p.is_file() and p.name not in ('bootstrap.py','README.md'):shutil.copy(p,out/p.name)
 for p in (ROOT/'core').iterdir():
  if p.suffix=='.py' or p.name in ['root20.json','anchor_angles20.json','anchor_isolation_certificate.json']:shutil.copy(p,out/p.name)
 for fn in ['cases20_meta.json','candidate20_map.json']:shutil.copy(ROOT/'enumeration'/fn,out/fn)
 for B in range(12,20):
  raw=gzip.decompress((ROOT/'enumeration'/f'survivors20_B{B}.bin.gz').read_bytes());tag,N,bb,n=struct.unpack_from('<4sIII',raw)
  if (tag,N,bb)!=(b'D20T',20,B):raise ValueError('Invalid packed input')
  rows=np.frombuffer(raw,dtype='<u4',offset=16).reshape((n,38-B))
  np.save(out/f'masks20_B{B}.npy',rows,allow_pickle=False)
 print('Created separate proposal workspace:',out)
 print('No accepting source or certificate in the proof directory was changed.')
if __name__=='__main__':main(sys.argv[1])
