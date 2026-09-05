"""Canonical geometry indices for the 119-variable 13-disk system."""
import json
from pathlib import Path
N=119
QSTART=0
FSTART=18
TID=62
WSTART=63
MUSTART=110
FIXED=0
QIDS=list(range(1,10))
FREEIDS=list(range(10,32))
EDGES=[tuple(e) for e in json.loads((Path(__file__).parent/'cover13_kkt_root_119d.json').read_text())['edges']]
