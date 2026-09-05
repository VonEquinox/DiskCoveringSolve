"""Independent standard-library replay of every remaining continuous domain."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,gzip,time,hashlib
import exact_force15 as G
import system15 as S
from candidate_global15 import local
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration'

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def traverse(B,faces,nodes,allow_local):
 assert isinstance(nodes,list) and nodes
 M=G.bounds(B,faces);ed=G.edges(B,faces);lo,hi=G.rootbox(M);seen=set();stack=[(0,lo,hi,0)];counts=Counter();mindual=minlocal=None;maxdepth=0
 while stack:
  idx,lo,hi,depth=stack.pop();assert type(idx)is int and 0<=idx<len(nodes) and idx not in seen;seen.add(idx);maxdepth=max(maxdepth,depth)
  node=nodes[idx];assert isinstance(node,dict);kind=node['kind'];box=G.tighten(M,lo,hi)
  if kind=='EMPTY':assert box is None;counts[kind]+=1;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   j=node['axis'];assert isinstance(node['mid'],str);cut=F(node['mid']);ch=node['children'];assert type(j)is int and 0<=j<B-1 and lo[j]<cut<hi[j]
   assert isinstance(ch,list) and len(ch)==2 and ch[0]!=ch[1] and all(type(a)is int for a in ch)
   lh=list(hi);lh[j]=cut;rl=list(lo);rl[j]=cut
   stack.extend([(ch[1],rl,list(hi),depth+1),(ch[0],list(lo),lh,depth+1)])
  elif kind=='DUAL':
   mar=G.check(B,ed,lo,hi,node);assert mar is not None and mar>0,(idx,mar);mindual=mar if mindual is None else min(mindual,mar)
  elif kind=='LOCAL':
   assert allow_local and B==S.B and [list(f) for f in faces]==[list(f) for f in S.FACES]
   mar=local(lo,hi);assert mar is not None and mar>0;minlocal=mar if minlocal is None else min(minlocal,mar)
  else:raise AssertionError(('unknown node/leaf',kind))
  counts[kind]+=1
 assert seen==set(range(len(nodes)))
 return counts,mindual,minlocal,maxdepth

def check_mapping(rows):
 cm=json.load(open(ENUM/'candidate15_map.json'));k=cm['residual_index'];assert type(k)is int and 0<=k<len(rows) and rows[k]==cm['row'] and rows[k]['B']==S.B
 p=cm['map'];assert len(p)==S.N and all(type(v)is int for v in p) and sorted(p)==list(range(S.N))
 assert sorted(p[S.B:])==list(range(S.B,S.N))
 assert any(all(p[i]==(sh+sgn*i)%S.B for i in range(S.B)) for sh in range(S.B) for sgn in [-1,1])
 assert sorted(tuple(sorted(p[i] for i in f)) for f in S.FACES)==sorted(tuple(f) for f in rows[k]['faces'])
 face_ix={tuple(f):i for i,f in enumerate(S.FACES)}
 mp=list(range(S.B+S.N))+[S.B+S.N+face_ix[tuple(f)] for f in S.ACTIVE]
 full={frozenset(e) for e in G.edges(S.B,S.FACES)}
 assert len(set(mp))==len(mp) and all(frozenset((mp[u],mp[v])) in full for u,v in S.EDGES)
 return k

def verify_candidate():
 cp=R/'candidate15_tree.json.gz'
 with gzip.open(cp,'rt') as f:c=json.load(f)
 assert c['version']==1 and c['B']==S.B and c['I']==S.N-S.B and c['faces']==[list(f) for f in S.FACES]
 assert c['bounds_matrix']==G.bounds(S.B,S.FACES)
 co,md,ml,dep=traverse(S.B,S.FACES,c['nodes'],True);assert co['LOCAL']>0
 out={'verified':True,'nodes':len(c['nodes']),'counts':dict(co),'minimum_force_margin':str(md),'minimum_local_margin':str(ml),'maxdepth':dep,'certificate_sha256':sha(cp)}
 print('CANDIDATE EXACTLY VERIFIED',dict(co),'force margin',float(md),'local margin',float(ml),flush=True)
 return out

def check_complete_indices(seen,rows,candidate,B=None):
 expected={k for k,row in enumerate(rows) if k!=candidate and (B is None or row['B']==B)}
 assert seen==expected,'incomplete/extra residual record set'

def verify():
 st=time.time();rows=json.load(open(ENUM/'metric15_residuals.json'));cand=check_mapping(rows);outcand=verify_candidate();rhash=sha(ENUM/'metric15_residuals.json')
 seen=set();counts=Counter();minmar=None;maxdepth=0;totalnodes=0;family={}
 for B in range(10,15):
  np=ENUM/f'noncandidate15_B{B}.jsonl.gz';fnodes=0;fcounts=Counter();fseen=set()
  with gzip.open(np,'rt') as f:
   h=json.loads(next(f));assert h['version']==1 and h['B']==B and h['I']==15-B and h['residual_sha256']==rhash and h['radius_upper']==str(G.RU)
   for line in f:
    case=json.loads(line);k=case['residual_index'];assert type(k)is int and 0<=k<len(rows) and k!=cand and k not in seen;seen.add(k);fseen.add(k);row=rows[k]
    assert row['B']==B and all(case[key]==row[key] for key in ['B','I','idx','faces','signature'])
    cs,m,l,d=traverse(B,row['faces'],case['nodes'],False);assert l is None
    counts.update(cs);fcounts.update(cs);minmar=m if minmar is None else (minmar if m is None else min(minmar,m));maxdepth=max(d,maxdepth);totalnodes+=len(case['nodes']);fnodes+=len(case['nodes'])
    if len(seen)%1000==0:print('NONCANDIDATE REPLAY',len(seen),'nodes',totalnodes,'seconds',time.time()-st,flush=True)
  check_complete_indices(fseen,rows,cand,B)
  family[str(B)]={'cases':len(fseen),'nodes':fnodes,'counts':dict(fcounts),'certificate_sha256':sha(np)}
 check_complete_indices(seen,rows,cand)
 out={'verified':True,'candidate':outcand,'noncandidate':{'cases':len(seen),'nodes':totalnodes,'counts':dict(counts),'minimum_force_margin':str(minmar),'minimum_force_margin_float':float(minmar),'maxdepth':maxdepth,'families':family},'root_sha256':sha(R/'root15.json'),'local_certificate_sha256':sha(R/'anchor_isolation_certificate.json'),'seconds':time.time()-st}
 (R/'forest15_verified.json').write_text(json.dumps(out,indent=2));print('ALL ANGLE TREES EXACTLY VERIFIED',len(seen),'cases',totalnodes,'nodes',flush=True);return out
if __name__=='__main__':verify()
