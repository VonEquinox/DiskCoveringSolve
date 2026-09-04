#!/usr/bin/env python3
"""Generate and exactly verify graph-energy certificates for all noncandidate
residual triangulations in the n=12 proof.

Floating point is used only to choose stress rows and the finite split tree.
Verification recomputes rationalized Kron coefficients and every leaf lower
bound using Fraction arithmetic and certified trigonometric minorants.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json,pickle,math,time,sys,heapq,functools
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from exact_core import PI_L,PI_U,TWO_PI_L,TWO_PI_U,exact_kron,tighten_sum_box,pair_intervals,linear_min
from fixed_trig import line_minorant
ROOT=Path(__file__).resolve().parent;DELTA=F(73891,100000);TARGET=F(13039536,10**8);QW=10**15;QD=10**18
CAND=(9,3,21015);HARD={(11,1,3101),(11,1,3128),(11,1,3161),(11,1,3216)}

def rationalize(w):
 w=np.maximum(np.asarray(w,float),0);w/=w.sum();raw=w*QW;nums=np.floor(raw).astype(object);nums=[int(x) for x in nums];rem=QW-sum(nums);frac=raw-np.floor(raw)
 for k in np.argsort(-frac)[:rem]:nums[int(k)]+=1
 assert min(nums)>=0 and sum(nums)==QW;return nums

def topology_data():
 metric=json.load(open(ROOT/'metric_exact_cert_safe.json'));disc=json.load(open(ROOT/'stress_discovery_safe.json'));dm={(z['B'],z['I'],z['idx']):z for z in disc};out=[]
 for key,b in metric['cases'].items():
  B,I=map(int,key.split(','));reps=pickle.load(open(ROOT/f'triangulations_B{B}_I{I}_orbits.pkl','rb'))
  for idx in b['residual_indices']:
   tup=(B,I,int(idx))
   if tup==CAND:continue
   z=dm[tup];weights=[z['edge_weights']]
   if tup in HARD:weights=json.load(open(ROOT/f'envelope_{B}_{I}_{idx}.json'))['weights']
   out.append({'B':B,'I':I,'idx':int(idx),'faces':[list(f) for f in reps[int(idx)]],'weights':weights})
 return out

def make_rows(z):
 B,I=z['B'],z['I'];faces=[tuple(f) for f in z['faces']];pairs=[(i,j) for i in range(B) for j in range(i+1,B)];recs=[]
 for w in z['weights']:
  nums=rationalize(w);D=exact_kron(B,I,faces,nums,QW);dn=[]
  for q in D:dn.append((q*QD).numerator//(q*QD).denominator)
  recs.append({'wnum':nums,'dnum':dn})
 return pairs,recs

def tighten_f(lo,hi):
 lo=list(lo);hi=list(hi);L=2*math.pi
 for _ in range(40):
  sl=sum(lo);sh=sum(hi)
  if sh<L-1e-13 or sl>L+1e-13:return None
  ch=False
  for i in range(len(lo)):
   nl=L-(sh-hi[i]);nh=L-(sl-lo[i])
   if nl>lo[i]:lo[i]=nl;ch=True
   if nh<hi[i]:hi[i]=nh;ch=True
   if lo[i]>hi[i]+1e-13:return None
  if not ch:break
 return lo,hi

def mincos(l,u):
 if l<=math.pi<=u:return -1.
 return min(math.cos(l),math.cos(u))
def line_f(l,u):
 a=(l+u)/2;h=(u-l)/2;A=math.sin(a);B=1-math.cos(a)-A*a+min(0.,mincos(l,u))*h*h/2-2e-13;return A,B

def linmin_f(c,b,lo,hi):
 x=np.where(c>=0,lo,hi);v=b+float(c@x);S=float(x.sum());L=2*math.pi
 if S<L:
  need=L-S
  for i in np.argsort(c):q=min(hi[i]-x[i],need);v+=c[i]*q;need-=q
 elif S>L:
  need=S-L
  for i in np.argsort(-c):q=min(x[i]-lo[i],need);v-=c[i]*q;need-=q
 return v

def bound_float(lo,hi,pairs,rows):
 z=tighten_f(lo,hi)
 if z is None:return 1e300,-1
 lo=np.array(z[0]);hi=np.array(z[1]);sl=lo.sum();sh=hi.sum();lines=[]
 for i,j in pairs:
  l0=lo[i:j].sum();u0=hi[i:j].sum();l=max(l0,2*math.pi-(sh-u0));u=min(u0,2*math.pi-(sl-l0));lines.append(line_f(l,u))
 best=-1e300;bi=-1
 for ri,row in enumerate(rows):
  c=np.zeros(len(lo));bb=0
  for n,(i,j),(a,b) in zip(row,pairs,lines):
   q=n/QD;bb+=q*b;c[i:j]+=q*a
  v=linmin_f(c,bb,lo,hi)
  if v>best:best=v;bi=ri
 return best,bi

def bound_exact(lo,hi,pairs,row):
 z=tighten_sum_box(lo,hi)
 if z is None:return None
 lo,hi=z;lines=[line_minorant(l,u,PI_L,PI_U) for l,u in pair_intervals(lo,hi,pairs)];c=[F(0)]*len(lo);b=F(0)
 for n,(i,j),(a,b0) in zip(row,pairs,lines):
  q=F(int(n),QD);b+=q*b0
  for k in range(i,j):c[k]+=q*a
 return b+linear_min(c,lo,hi)

def generate_one(z,maxnodes=2000000,margin=2e-9):
 pairs,recs=make_rows(z);rows=[r['dnum'] for r in recs];B=z['B'];nodes=[];boxes=[];stack=[((0.,)*B,(float(DELTA),)*B,None,None)];t=time.time()
 while stack:
  lo,hi,par,side=stack.pop();idx=len(nodes);nodes.append(None);boxes.append((lo,hi))
  if par is not None:nodes[par]['child'][side]=idx
  v,ri=bound_float(np.array(lo),np.array(hi),pairs,rows)
  if v>1e200:nodes[idx]={'empty':1};continue
  if v>=float(TARGET)+margin:nodes[idx]={'leaf':ri};continue
  widths=[hi[i]-lo[i] for i in range(B)];k=max(range(B),key=lambda i:(widths[i],-i));m=(lo[k]+hi[k])/2
  nodes[idx]={'split':k,'child':[None,None]};lhi=list(hi);lhi[k]=m;rlo=list(lo);rlo[k]=m
  stack.append((tuple(rlo),hi,idx,1));stack.append((lo,tuple(lhi),idx,0))
  if len(nodes)>maxnodes:raise RuntimeError(('too many',z['B'],z['I'],z['idx'],len(nodes),v))
 cert={k:z[k] for k in ['B','I','idx','faces']};cert.update({'qweight':QW,'qcoef':QD,'records':recs,'nodes':nodes})
 print('GENONE',(z['B'],z['I'],z['idx']),'records',len(recs),'nodes',len(nodes),'leaves',sum('leaf'in n for n in nodes),'sec',time.time()-t,flush=True);return cert

def reconstruct(nodes,B):
 boxes=[None]*len(nodes);boxes[0]=((F(0),)*B,(DELTA,)*B);leaves=[]
 for i,n in enumerate(nodes):
  assert boxes[i] is not None;lo,hi=boxes[i]
  if 'split'in n:
   k=int(n['split']);m=(lo[k]+hi[k])/2
   for side,j in enumerate(n['child']):
    l=list(lo);h=list(hi)
    if side==0:h[k]=m
    else:l[k]=m
    assert boxes[int(j)] is None;boxes[int(j)]=(tuple(l),tuple(h))
  elif 'leaf'in n:leaves.append(i)
  elif 'empty'in n:pass
  else:raise AssertionError(n)
 return boxes,leaves

def verify_one(z):
 B,I=z['B'],z['I'];faces=[tuple(f) for f in z['faces']];pairs=[(i,j) for i in range(B) for j in range(i+1,B)];assert z['qweight']==QW and z['qcoef']==QD
 rows=[]
 for rec in z['records']:
  nums=list(map(int,rec['wnum']));assert min(nums)>=0 and sum(nums)==QW;D=exact_kron(B,I,faces,nums,QW);stored=list(map(int,rec['dnum']));assert len(stored)==len(D)
  for n,q in zip(stored,D):assert F(n,QD)<=q
  rows.append(stored)
 boxes,leaves=reconstruct(z['nodes'],B);mm=None
 for i in leaves:
  ri=int(z['nodes'][i]['leaf']);v=bound_exact(*boxes[i],pairs,rows[ri]);assert v is not None and v>=TARGET,(B,I,z['idx'],i,float(v-TARGET));mar=v-TARGET;mm=mar if mm is None or mar<mm else mm
 for i,n in enumerate(z['nodes']):
  if 'empty'in n:assert tighten_sum_box(*boxes[i]) is None
 return {'B':B,'I':I,'idx':z['idx'],'nodes':len(z['nodes']),'leaves':len(leaves),'records':len(rows),'min_margin':float(mm)}

def generate_all(path,only=None):
 arr=[]
 for z in topology_data():
  if only and (z['B'],z['I'],z['idx'])!=only:continue
  arr.append(generate_one(z))
 json.dump({'version':1,'target':[TARGET.numerator,TARGET.denominator],'delta':[DELTA.numerator,DELTA.denominator],'topologies':arr},open(path,'w'),separators=(',',':'));print('saved',path)
def verify_all(path):
 D=json.load(open(path));assert F(*D['target'])==TARGET and F(*D['delta'])==DELTA;out=[];t=time.time()
 # Independently identify the expected 131 noncandidate residual orbits and their faces.
 metric=json.load(open(ROOT/'metric_exact_cert_safe.json'));expected={}
 for key,b in metric['cases'].items():
  B,I=map(int,key.split(','));reps=pickle.load(open(ROOT/f'triangulations_B{B}_I{I}_orbits.pkl','rb'))
  for idx in b['residual_indices']:
   tup=(B,I,int(idx))
   if tup!=CAND:expected[tup]=[list(f) for f in reps[int(idx)]]
 actual={(int(z['B']),int(z['I']),int(z['idx'])):z['faces'] for z in D['topologies']}
 assert len(actual)==len(D['topologies']) and actual==expected and len(actual)==131
 for z in D['topologies']:
  q=verify_one(z);out.append(q);print('VERONE',q,flush=True)
 rep={'verified':True,'partition_verified':True,'topologies':len(out),'nodes':sum(x['nodes'] for x in out),'leaves':sum(x['leaves'] for x in out),'min_margin':min(x['min_margin'] for x in out),'seconds':time.time()-t,'details':out};json.dump(rep,open(ROOT/'residual_energy_safe_verified.json','w'),indent=2);print(json.dumps({k:v for k,v in rep.items() if k!='details'},indent=2))
if __name__=='__main__':
 mode=sys.argv[1];path=sys.argv[2]
 only=tuple(map(int,sys.argv[3:6])) if len(sys.argv)>=6 else None
 generate_all(path,only) if mode=='generate' else verify_all(path)
