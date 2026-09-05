"""Exact partition linkage; C++ audits every state and every Farkas inequality."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
from math import comb
import json,hashlib,time,subprocess
from exact_arcs import sincos_turn,SCALE
import exact_force15 as G
import system15 as S
if not __debug__:raise RuntimeError('Run without -O')
R=Path(__file__).resolve().parent;ENUM=R.parent/'enumeration';WORK=R.parent/'replay_workspace';REPORTS=R.parent/'reports'
@lru_cache(maxsize=None)
def labelled(B,I):
 assert B>=2 and I>=0
 if B==2:return int(I==0)
 total=sum(comb(I,j)*labelled(k,j)*labelled(B-k+1,I-j) for k in range(2,B) for j in range(I+1))
 if I:total+=I*(labelled(B+1,I-1)-sum(comb(I-1,j)*labelled(3,j)*labelled(B,I-1-j) for j in range(I)))
 return total

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def check_constants():
 d,z,Y,rho=S.load();t=z[S.TID];rho=F(1,10**80);A,C,D,RU=G.A,G.C,G.D,G.RU
 assert 0<t-rho<t+rho<RU**2<F(1,4)
 _,s1=sincos_turn(F(A,2*D));_,s2=sincos_turn(F(C,2*D));_,s9=sincos_turn(F(1,18))
 assert RU<F(s1[0],SCALE) and 2*RU<F(s2[0],SCALE) and RU<F(s9[0],SCALE)
 assert 0<2*A<C<D//2
 return {'radius_upper':str(RU),'gap_cap_turns':str(F(A,D)),'four_rod_cap_turns':str(F(C,D)),'minimum_boundary_cells':10,'sine_margins':[str(F(s1[0],SCALE)-RU),str(F(s2[0],SCALE)-2*RU),str(F(s9[0],SCALE)-RU)]}

def verify():
 st=time.time();constants=check_constants();cert=json.load(open(ENUM/'metric15_certificate.json'));triples=list(combinations(range(15),3));allres=[];summary={};tot=0;mass={}
 for B in range(10,15):
  fam=f'B{B}_I{15-B}';p=WORK/(fam+'.txt');c=cert[fam];assert sha(p)==c['source_txt_sha256']
  cp=WORK/(fam+'.classes.txt');rp=WORK/(fam+'.metric_residual.txt');sp=WORK/(fam+'.metric_summary.txt')
  cache={tuple(map(int,k.split(','))) if k else ():v for k,v in c['certified_signatures'].items()};survive={tuple(map(int,k.split(','))) if k else () for k in c['surviving_signatures']}
  assert not(set(cache)&survive)
  with open(cp,'w') as f:
   print(B,len(cache)+len(survive),file=f)
   for key in sorted(set(cache)|survive):
    assert list(key)==sorted(set(key));w=cache.get(key);line=[int(w is not None),len(key),*key]
    if w is not None:
     assert len(w)==B+len(key) and all(type(v)is int and v>=0 for v in w);line+=w
    print(*line,file=f)
  rp.unlink(missing_ok=True);sp.unlink(missing_ok=True)
  subprocess.run([str(WORK/'metric_audit15'),str(p),str(cp),str(rp),str(sp)],check=True)
  b,i,n,dr,fa,re=map(int,sp.read_text().split());assert (b,i)==(B,15-B);counts={'direct':dr,'farkas':fa,'residual':re};assert counts==c['counts'];summary[fam]=counts;tot+=n;indices=[]
  with open(rp) as f:
   for line in f:
    ar=list(map(int,line.split()));idx,k=ar[:2];key=ar[2:2+k];ids=ar[2+k:];assert len(ids)==28-B and key==sorted(set(key)) and idx>=0 and (not indices or idx>indices[-1]);indices.append(idx)
    row={'B':B,'I':15-B,'idx':idx,'faces':[list(triples[j]) for j in ids],'signature':key}
    assert G.signature(B,row['faces'])==tuple(key)
    allres.append(row)
  assert len(indices)==re and indices==c['residual_indices']
  audit=json.load(open(REPORTS/f'audit{B}.json'));expected=labelled(B,15-B)
  assert audit['verified'] and audit['B']==B and audit['I']==15-B and audit['representatives']==n and audit['orbit_mass']==expected
  mass[fam]=expected
  print('PARTITION FAMILY LINKED',fam,counts,'seconds',time.time()-st,flush=True)
 assert set(cert)==set(summary) and allres==json.load(open(ENUM/'metric15_residuals.json'))
 out={'verified':True,'constants':constants,'families':summary,'total':tot,'excluded':sum(z['direct']+z['farkas'] for z in summary.values()),'residual':len(allres),'labelled_totals_from_recurrence':mass,'seconds':time.time()-st}
 (R/'partition15_verified.json').write_text(json.dumps(out,indent=2));print('PARTITION VERIFIED',tot,'residual',len(allres),flush=True);return out
if __name__=='__main__':verify()
