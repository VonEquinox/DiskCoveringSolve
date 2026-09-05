"""Exact cross-stage graph, radius, and angular constant interfaces."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
import system20 as S
import interval20 as I
import domain20 as D
import compact_force20 as C
from exact_arcs import sincos_turn,SCALE
from cases20 import Store,ENUM
if not __debug__:raise RuntimeError('Run without -O')
def verify():
 assert S.N==D.N==C.N==20 and S.B==13 and D.V==C.V==58
 assert (S.G,S.P,S.TID,S.M,S.C,S.D)==(102,103,102,83,95,198)
 assert len(S.ACTIVE)==19 and len(S.FACES)==25 and set(S.ACTIVE)<=set(S.FACES)
 assert len(set(S.ACTIVE))==19 and len(set(S.FACES))==25
 ed=[(i,S.B+i) for i in range(S.B)]+[(i,S.B+(i+1)%S.B) for i in range(S.B)]+[(S.B+S.N+k,S.B+c) for k,fa in enumerate(S.ACTIVE) for c in fa]
 assert ed==S.EDGES and len(S.CONS)==S.C
 dat,z,iv,p=I.rootbox();assert dat['N']==198 and dat['all_faces']==[list(f) for f in S.FACES]
 assert C.RU==D.RU and 0<iv[S.TID][0]<iv[S.TID][1]<D.RU**2<F(1,4)
 caps=[F(D.A,D.D),F(D.C4,D.D),F(D.C6,D.D)];assert (D.A,D.C4,D.C6,D.D)==(873043,1822179,3019184,10000000)
 assert 0<caps[0]<caps[1]<caps[2]<F(1,2) and 11*caps[0]<1
 margins=[]
 for k,cap in enumerate(caps,1):
  _,si=sincos_turn(cap/2);mar=F(si[0],SCALE)-k*D.RU;assert mar>0;margins.append(str(mar))
 store=Store();meta=json.load(open(ENUM/'cases20_meta.json'));assert meta['N']==20
 assert [(r['B'],r['I'],r['orbits']) for r in meta['families']]==[(b,20-b,store.counts[b]) for b in range(12,20)]
 c=json.load(open(ENUM/'candidate20_map.json'));B=c['B'];idx=c['idx'];mp=c['map'];fs=store.faces(B,idx)
 assert B==S.B and type(idx)is int and len(mp)==20 and sorted(mp)==list(range(20))
 assert all(mp[i]==(mp[0]+i)%B for i in range(B)) and set(mp[B:])==set(range(B,20))
 assert c['faces']==fs and {tuple(sorted(mp[v] for v in f)) for f in S.FACES}==set(map(tuple,fs))
 lookup={tuple(f):j for j,f in enumerate(fs)}
 node_map=[mp[i] for i in range(B)]+[B+mp[i] for i in range(20)]+[B+20+lookup[tuple(sorted(mp[v] for v in f))] for f in S.ACTIVE]
 assert len(node_map)==52 and len(set(node_map))==52
 target={frozenset(e) for e in D.edges(B,fs)}
 assert all(frozenset((node_map[u],node_map[v])) in target for u,v in S.EDGES)
 out={'verified':True,'radius_upper':str(D.RU),'turn_caps':list(map(str,caps)),'sine_margins':margins,'minimum_boundary_cells':12,'families':[[b,20-b] for b in range(12,20)],'candidate_B':B,'candidate_index':idx,'center_map':mp,'active_node_embedding':node_map,'root_sha256':hashlib.sha256((S.R/'root20.json').read_bytes()).hexdigest()}
 (S.R/'interfaces20_verified.json').write_text(json.dumps(out,indent=2));print('EXACT INTERFACES VERIFIED; candidate',B,idx,flush=True);return out
if __name__=='__main__':verify()
