"""Replay every noncandidate angle tree using only exact standard-library arithmetic."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,gzip,time,hashlib
import exact_force19 as G
import geometry19 as S
if not __debug__:raise RuntimeError('Assertions must be enabled')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def traverse(B,faces,nodes):
 M=G.bounds(B,faces);ed=G.edges(B,faces);lo,hi=G.rootbox(M);seen=set();stack=[(0,lo,hi,0)];counts=Counter();minimum=None;maxdepth=0
 while stack:
  idx,lo,hi,depth=stack.pop();assert type(idx)is int and 0<=idx<len(nodes) and idx not in seen;seen.add(idx);maxdepth=max(maxdepth,depth)
  node=nodes[idx];kind=node['kind'];box=G.tighten(M,lo,hi)
  if kind=='EMPTY':assert box is None;counts[kind]+=1;continue
  assert box is not None;lo,hi=box
  if kind=='SPLIT':
   j=node['axis'];cut=F(node['mid']);ch=node['children'];assert type(j)is int and 0<=j<B-1 and lo[j]<cut<hi[j]
   assert type(ch)is list and len(ch)==2 and ch[0]!=ch[1] and all(type(a)is int for a in ch)
   lh=list(hi);lh[j]=cut;rl=list(lo);rl[j]=cut
   stack.extend([(ch[1],rl,list(hi),depth+1),(ch[0],list(lo),lh,depth+1)])
  elif kind=='DUAL':
   mar=G.check(B,ed,lo,hi,node);assert mar is not None and mar>0,(idx,mar);minimum=mar if minimum is None else min(minimum,mar)
  else:raise AssertionError(('unknown or impermissible leaf',kind))
  counts[kind]+=1
 assert seen==set(range(len(nodes)))
 return counts,minimum,maxdepth

def check_mapping(rows,cm=None):
 if cm is None:cm=json.load(open(ENUM/'candidate19_map.json'))
 k=cm['residual_index'];assert type(k)is int and 0<=k<len(rows) and rows[k]==cm['row'] and rows[k]['B']==S.B
 p=cm['map'];assert len(p)==S.N and all(type(v)is int for v in p) and sorted(p)==list(range(S.N))
 assert sorted(p[S.B:])==list(range(S.B,S.N))
 assert any(all(p[i]==(sh+sgn*i)%S.B for i in range(S.B)) for sh in range(S.B) for sgn in [-1,1])
 assert sorted(tuple(sorted(p[i] for i in f)) for f in S.FACES)==sorted(tuple(f) for f in rows[k]['faces'])
 full={frozenset(e) for e in G.edges(S.B,S.FACES)}
 assert all(frozenset(e) in full for e in S.EDGES)
 return k

def complete_case_set(seen,n,candidate):
 assert candidate not in seen and seen==set(range(n))-{candidate},'missing or extra topology'

def verify():
 (R/'forest19_verified.json').unlink(missing_ok=True);st=time.time();rows=json.load(open(ENUM/'metric19_residuals.json'));cand=check_mapping(rows)
 path=ENUM/'noncandidate19_trees.jsonl.gz';seen=set();counts=Counter();minimum=None;maxdepth=0;totalnodes=0;families={}
 with gzip.open(path,'rt',encoding='utf8') as f:
  head=json.loads(f.readline());assert head['version']==1 and head['residual_sha256']==sha(ENUM/'metric19_residuals.json') and head['radius_upper']==str(G.RU) and head['cases']==len(rows)-1
  for line in f:
   case=json.loads(line);k=case['residual_index'];assert type(k)is int and 0<=k<len(rows) and k!=cand and k not in seen;seen.add(k);row=rows[k]
   assert all(case[key]==row[key] for key in ['B','I','idx','faces','signature','stabilizer'])
   co,mar,depth=traverse(row['B'],row['faces'],case['nodes']);assert mar is not None and mar>0
   counts.update(co);minimum=mar if minimum is None else min(minimum,mar);maxdepth=max(maxdepth,depth);totalnodes+=len(case['nodes'])
   b=str(row['B']);fa=families.setdefault(b,{'cases':0,'nodes':0,'DUAL':0,'SPLIT':0,'EMPTY':0});fa['cases']+=1;fa['nodes']+=len(case['nodes'])
   for key in ['DUAL','SPLIT','EMPTY']:fa[key]+=co[key]
   if len(seen)%5000==0:print('EXACT NONCANDIDATE REPLAY',len(seen),'nodes',totalnodes,'seconds',time.time()-st,flush=True)
 complete_case_set(seen,len(rows),cand)
 assert minimum>0 and counts['DUAL']+counts['SPLIT']+counts['EMPTY']==totalnodes
 out={'verified':True,'candidate_index':cand,'candidate_cases':1,'candidate_method':'global rational Dirichlet-energy strict convexity, no local leaves','noncandidate_cases':len(seen),'nodes':totalnodes,'counts':dict(counts),'minimum_force_margin':str(minimum),'minimum_force_margin_float':float(minimum),'maxdepth':maxdepth,'families':families,'unresolved_leaves':0,'certificate_sha256':sha(path),'residual_sha256':sha(ENUM/'metric19_residuals.json'),'seconds':time.time()-st}
 (R/'forest19_verified.json').write_text(json.dumps(out,indent=2));print('ALL NONCANDIDATE ANGLE DOMAINS VERIFIED',len(seen),dict(counts),flush=True);return out
if __name__=='__main__':verify()
