"""Independent standard-library replay of the complete remaining angle space."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,gzip,time,hashlib
import exact_force14 as G
import system14 as S
from candidate_global14 import local
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def traverse(B,faces,nodes,allow_local):
 M=G.bounds(B,faces);ed=G.edges(B,faces);lo,hi=G.rootbox(M);seen=set();stack=[(0,lo,hi,0)];counts=Counter();mindual=minlocal=None;maxdepth=0
 while stack:
  idx,lo,hi,depth=stack.pop();assert type(idx)is int and 0<=idx<len(nodes) and idx not in seen;seen.add(idx);maxdepth=max(maxdepth,depth)
  node=nodes[idx];kind=node['kind'];box=G.tighten(M,lo,hi)
  if kind=='EMPTY':assert box is None;counts[kind]+=1;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   j=node['axis'];cut=F(node['mid']);ch=node['children'];assert type(j)is int and 0<=j<B-1 and lo[j]<cut<hi[j]
   assert isinstance(ch,list) and len(ch)==2 and ch[0]!=ch[1] and all(type(a)is int for a in ch)
   lh=list(hi);lh[j]=cut;rl=list(lo);rl[j]=cut
   stack.extend([(ch[1],rl,list(hi),depth+1),(ch[0],list(lo),lh,depth+1)])
  elif kind=='DUAL':
   mar=G.check(B,ed,lo,hi,node);assert mar is not None and mar>0,(idx,mar);mindual=mar if mindual is None else min(mindual,mar)
  elif kind=='LOCAL':
   assert allow_local and B==10 and [list(f) for f in faces]==[list(f) for f in S.FACES]
   mar=local(lo,hi);assert mar is not None and mar>0;minlocal=mar if minlocal is None else min(minlocal,mar)
  else:raise AssertionError(('unknown leaf',kind))
  counts[kind]+=1
 assert seen==set(range(len(nodes)))
 return counts,mindual,minlocal,maxdepth

def check_mapping(rows):
 cm=json.load(open(ENUM/'candidate14_map.json'));k=cm['residual_index'];assert type(k)is int and 0<=k<len(rows) and rows[k]==cm['row'] and rows[k]['B']==10
 p=cm['map'];assert len(p)==14 and sorted(p)==list(range(14))
 assert sorted(p[10:])==[10,11,12,13]
 assert any(all(p[i]==(sh+sgn*i)%10 for i in range(10)) for sh in range(10) for sgn in [-1,1])
 assert sorted(tuple(sorted(p[i] for i in f)) for f in S.FACES)==sorted(tuple(f) for f in rows[k]['faces'])
 # Exact embedding of the positive-weight candidate graph into the full graph.
 face_ix={tuple(f):i for i,f in enumerate(S.FACES)}
 mp=list(range(24))+[24+face_ix[tuple(f)] for f in S.ACTIVE]
 full={frozenset(e) for e in G.edges(10,S.FACES)}
 assert all(frozenset((mp[u],mp[v])) in full for u,v in S.EDGES)
 return k

def verify():
 st=time.time();rows=json.load(open(ENUM/'metric14_residuals.json'));cand=check_mapping(rows)
 cp=R/'candidate14_tree.json.gz'
 with gzip.open(cp,'rt') as f:c=json.load(f)
 assert c['version']==1 and c['faces']==[list(f) for f in S.FACES]
 co,md,ml,dep=traverse(10,S.FACES,c['nodes'],True)
 assert co['LOCAL']>0
 outcand={'verified':True,'nodes':len(c['nodes']),'counts':dict(co),'minimum_force_margin':str(md),'minimum_local_margin':str(ml),'maxdepth':dep,'certificate_sha256':sha(cp)}
 print('CANDIDATE EXACTLY VERIFIED',dict(co),'force margin',float(md),'local margin',float(ml),flush=True)
 np=ENUM/'noncandidate14_trees.json.gz'
 with gzip.open(np,'rt') as f:nc=json.load(f)
 assert nc['version']==1 and nc['residual_sha256']==sha(ENUM/'metric14_residuals.json') and nc['radius_upper']==str(G.RU)
 seen=set();counts=Counter();minmar=None;maxdepth=0;totalnodes=0
 for case in nc['cases']:
  k=case['residual_index'];assert type(k)is int and 0<=k<len(rows) and k!=cand and k not in seen;seen.add(k);row=rows[k]
  assert all(case[key]==row[key] for key in ['B','I','idx','faces','signature'])
  cs,m,l,d=traverse(row['B'],row['faces'],case['nodes'],False);assert l is None
  counts.update(cs);minmar=m if minmar is None else (minmar if m is None else min(minmar,m));maxdepth=max(d,maxdepth);totalnodes+=len(case['nodes'])
  if len(seen)%500==0:print('NONCANDIDATE REPLAY',len(seen),'nodes',totalnodes,'seconds',time.time()-st,flush=True)
 assert seen==set(range(len(rows)))-{cand} and len(seen)==7231
 out={'verified':True,'candidate':outcand,'noncandidate':{'cases':len(seen),'nodes':totalnodes,'counts':dict(counts),'minimum_force_margin':str(minmar),'minimum_force_margin_float':float(minmar),'maxdepth':maxdepth,'certificate_sha256':sha(np)},'root_sha256':sha(R/'root14.json'),'local_certificate_sha256':sha(R/'anchor_isolation_certificate.json'),'seconds':time.time()-st}
 (R/'forest14_verified.json').write_text(json.dumps(out,indent=2));print('ALL ANGLE TREES EXACTLY VERIFIED',out['noncandidate'],flush=True);return out
if __name__=='__main__':verify()
