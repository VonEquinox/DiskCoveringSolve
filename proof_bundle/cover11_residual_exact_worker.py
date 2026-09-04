#!/usr/bin/env python3
"""Generate and replay exact T+ certificates for residual noncandidate topologies.

Numerical adaptive stresses are proposals only.  Each original graph weight is
rounded to an exact dyadic simplex point.  The anchor conductances are then
recomputed by Fraction star-mesh elimination and rounded downward to a common
dyadic denominator.  Every leaf is independently checked with rational pi and
Taylor bounds from cover11_candidate_exact_core.
"""
from __future__ import annotations
import argparse, json, math, os, sys, time, traceback
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core

RALL=json.load(open('/mnt/data/cover11_residual_cert_T14438.json'))
M=json.load(open('/mnt/data/cover11_metric_cert_T14438.json'))
FM={}
for key,c in M['cases'].items():
    B,I=map(int,key.split(','))
    for z in c['residual']:
        FM[(B,I,int(z['idx']))]=[tuple(map(int,f)) for f in z['faces']]

QWBITS=100; QW=1<<QWBITS
QDBITS=70; QD=1<<QDBITS
LMIXBITS=52; LQ=1<<LMIXBITS
TARGET=core.TARGET


def graph(B:int,I:int,faces):
    edges=[]
    for i in range(B):
        edges += [(i,B+i),(i,B+(i+1)%B)]
    off=2*B+I
    for k,fa in enumerate(faces):
        edges += [(off+k,B+v) for v in fa]
    n=off+len(faces)
    return edges,n


def round_simplex(w):
    w=np.maximum(np.asarray(w,float),0.0)
    w=w/w.sum();nums=[int(math.floor(float(x)*QW)) for x in w]
    p=int(np.argmax(w));nums[p]+=QW-sum(nums)
    assert min(nums)>=0 and sum(nums)==QW
    return nums


def exact_D(B,I,faces,nums):
    edges,n=graph(B,I,faces);assert len(nums)==len(edges)
    adj={i:{} for i in range(n)}
    def add(a,b,z):
        if a==b or z==0:return
        adj[a][b]=adj[a].get(b,F(0))+z
        adj[b][a]=adj[b].get(a,F(0))+z
    for (a,b),num in zip(edges,nums):add(a,b,F(num,QW))
    free=set(range(B,n))
    # Minimum-degree elimination is mathematically identical to a Schur
    # complement but avoids catastrophic Fraction fill-in on some graphs.
    while free:
        v=min(free,key=lambda z:(sum(1 for x in adj.get(z,{}).values() if x),z))
        free.remove(v)
        if v not in adj:continue
        nei=[(u,z) for u,z in adj[v].items() if z]
        s=sum((z for _,z in nei),F(0))
        if s:
            for i,(a,x) in enumerate(nei):
                for b,y in nei[i+1:]:add(a,b,x*y/s)
            for a,_ in nei:adj[a].pop(v,None)
        del adj[v]
    pairs=[(i,j) for i in range(B) for j in range(i+1,B)]
    exact=[2*adj.get(i,{}).get(j,F(0)) for i,j in pairs]
    dn=[]
    for z in exact:
        assert z>=0
        n0=(z*QD).numerator//(z*QD).denominator
        assert F(n0,QD)<=z
        dn.append(int(n0))
    return pairs,dn,exact


def exactify_stresses(B,I,faces,weights):
    rec=[];rows=[];mindown=None
    for w in weights:
        nums=round_simplex(w);pairs,dn,ex=exact_D(B,I,faces,nums)
        rows.append([F(n,QD) for n in dn])
        for n,z in zip(dn,ex):
            gap=z-F(n,QD)
            mindown=gap if mindown is None or gap<mindown else mindown
        rec.append({'weight_nums':[str(n) for n in nums],'Dnums':[str(n) for n in dn]})
    return pairs,rows,rec,mindown


def reconstruct_old(entry):
    B=int(entry['B']);nodes=entry['nodes'];boxes=[None]*len(nodes)
    boxes[0]=((F(0),)*B,(core.DELTA,)*B);leaves=[]
    for idx,n in enumerate(nodes):
        assert boxes[idx] is not None,(idx,'unreached')
        lo,hi=boxes[idx]
        if 'split' in n:
            k=int(n['split']);m=(lo[k]+hi[k])/2
            h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
            c0,c1=map(int,n['child'])
            assert boxes[c0] is None and boxes[c1] is None
            boxes[c0]=(lo,tuple(h0));boxes[c1]=(tuple(l1),hi)
        elif 'leaf' in n:leaves.append((idx,lo,hi))
        elif 'empty' in n:
            assert core.tighten_exact(lo,hi) is None
        else:raise AssertionError((idx,n))
    return leaves


def affine_data_exact(lo,hi,pairs,rows):
    th=core.tighten_exact(lo,hi)
    if th is None:return None
    lo,hi=th
    lines=[core.line_exact(l,u) for l,u in core.pair_intervals_exact(lo,hi,pairs)]
    A=[];bb=[]
    for row in rows:
        c=[F(0)]*len(lo);b=F(0)
        for coef,(i,j),(a,b0) in zip(row,pairs,lines):
            b+=coef*b0;da=coef*a
            for h in range(i,j):c[h]+=da
        A.append(c);bb.append(b)
    return lo,hi,A,bb


def row_exact_bound(c,b,lo,hi):
    return core.linear_min_exact(c,b,lo,hi)


def propose_mix(A,bb,lo,hi):
    K=len(A);B=len(lo)
    Af=np.array([[float(x) for x in r] for r in A],float)
    bf=np.array([float(x) for x in bb],float)
    lof=np.array([float(x) for x in lo]);hif=np.array([float(x) for x in hi])
    # Equality at central 2*pi is only used to propose the dual. Exact replay
    # handles the rigorous [S_L,S_U] enclosure.
    c=np.r_[np.zeros(B),1.0];scale=1e6
    Aub=np.column_stack([Af,-np.ones(K)])*scale;bub=-bf*scale
    rr=linprog(c,A_ub=Aub,b_ub=bub,
               A_eq=np.r_[np.ones(B),0.0][None,:],b_eq=[2*math.pi],
               bounds=[(lof[i],hif[i]) for i in range(B)]+[(None,None)],
               method='highs-ds',options={'primal_feasibility_tolerance':1e-10,
                 'dual_feasibility_tolerance':1e-10})
    if not rr.success:return None
    lam=np.maximum(-rr.ineqlin.marginals[:K]*scale,0.0)
    if lam.sum()<=0:return None
    lam/=lam.sum();nums=[int(math.floor(float(x)*LQ)) for x in lam]
    p=int(np.argmax(lam));nums[p]+=LQ-sum(nums)
    if min(nums)<0 or sum(nums)!=LQ:return None
    return nums


def mix_affine(A,bb,nums):
    c=[]
    for j in range(len(A[0])):
        c.append(sum((F(nums[k],LQ)*A[k][j] for k in range(len(A))),F(0)))
    b=sum((F(nums[k],LQ)*bb[k] for k in range(len(A))),F(0))
    return c,b


def witness(lo,hi,pairs,rows):
    dat=affine_data_exact(lo,hi,pairs,rows)
    if dat is None:return {'empty':1},None
    lo1,hi1,A,bb=dat
    vals=[row_exact_bound(A[k],bb[k],lo1,hi1) for k in range(len(rows))]
    k=max(range(len(vals)),key=lambda i:vals[i])
    if vals[k]>=TARGET:return {'kind':'row','row':k},vals[k]
    nums=propose_mix(A,bb,lo1,hi1)
    if nums is not None:
        c,b=mix_affine(A,bb,nums);v=row_exact_bound(c,b,lo1,hi1)
        if v>=TARGET:return {'kind':'mix','nums':[str(n) for n in nums],'den':str(LQ)},v
    return None,max(vals)


def split_index(lo,hi):
    B=len(lo);lev=[(i+1)*(B-i) for i in range(B)]
    return max(range(B),key=lambda i:((hi[i]-lo[i])*lev[i],-i))


def build_tree(lo0,hi0,pairs,rows,maxdepth=30,maxnodes=200000):
    nodes=[];stack=[(lo0,hi0,None,None,0)];mn=None
    while stack:
        lo,hi,parent,side,dep=stack.pop();idx=len(nodes);nodes.append(None)
        if parent is not None:nodes[parent]['child'][side]=idx
        wit,v=witness(lo,hi,pairs,rows)
        if wit is not None:
            if 'empty' not in wit:
                mar=v-TARGET;assert mar>=0
                mn=mar if mn is None or mar<mn else mn
            nodes[idx]=wit;continue
        if dep>=maxdepth:raise RuntimeError(('hard',dep,float(v-TARGET),[float(x) for x in lo],[float(x) for x in hi]))
        if len(nodes)>=maxnodes:raise RuntimeError(('maxnodes',len(nodes),float(v-TARGET)))
        k=split_index(lo,hi);m=(lo[k]+hi[k])/2
        nodes[idx]={'split':k,'child':[None,None]}
        h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
        stack.append((tuple(l1),hi,idx,1,dep+1));stack.append((lo,tuple(h0),idx,0,dep+1))
    return nodes,mn


def verify_tree(lo0,hi0,nodes,pairs,rows):
    boxes=[None]*len(nodes);boxes[0]=(lo0,hi0);mn=None;nleaf=nempty=0
    for idx,n in enumerate(nodes):
        assert boxes[idx] is not None,(idx,'unreached')
        lo,hi=boxes[idx]
        if 'split' in n:
            k=int(n['split']);m=(lo[k]+hi[k])/2
            h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
            c0,c1=map(int,n['child']);assert boxes[c0] is None and boxes[c1] is None
            boxes[c0]=(lo,tuple(h0));boxes[c1]=(tuple(l1),hi)
        elif 'empty' in n:
            assert core.tighten_exact(lo,hi) is None;nempty+=1
        elif n.get('kind') in ('row','mix'):
            dat=affine_data_exact(lo,hi,pairs,rows);assert dat is not None
            l,h,A,bb=dat
            if n['kind']=='row':c,b=A[int(n['row'])],bb[int(n['row'])]
            else:
                nums=list(map(int,n['nums']));den=int(n['den']);assert den==LQ and min(nums)>=0 and sum(nums)==den
                c,b=mix_affine(A,bb,nums)
            v=row_exact_bound(c,b,l,h);mar=v-TARGET
            assert mar>=0,(idx,n['kind'],float(mar))
            mn=mar if mn is None or mar<mn else mn;nleaf+=1
        else:raise AssertionError((idx,n))
    return nleaf,nempty,mn


def process_entry(k:int,verify=True):
    entry=RALL['entries'][k];B=int(entry['B']);I=int(entry['I']);idx=int(entry['idx'])
    af=Path(f'/mnt/data/resadapt_{k:02d}_{B}_{idx}.json')
    if not af.exists():raise FileNotFoundError(af)
    A=json.load(open(af));assert int(A['B'])==B and int(A['I'])==I and int(A['idx'])==idx
    faces=FM[(B,I,idx)];assert [tuple(f) for f in A['faces']]==faces
    pairs,rows,srecs,mindown=exactify_stresses(B,I,faces,A['weights'])
    roots=[];tot=0;mn=None;t=time.time();oldleaves=reconstruct_old(entry)
    for ri,(old,lo,hi) in enumerate(oldleaves):
        ns,mm=build_tree(lo,hi,pairs,rows);tot+=len(ns)
        mn=mm if mn is None or (mm is not None and mm<mn) else mn
        roots.append({'old_node':old,'nodes':ns})
    cert={'version':1,'entry_index':k,'B':B,'I':I,'idx':idx,
          'faces':[list(f) for f in faces],'pairs':[list(p) for p in pairs],
          'qwb':QWBITS,'qdb':QDBITS,'lmixb':LMIXBITS,'stresses':srecs,
          'roots':roots,'nodes':tot,
          'min_downward_gap':[mindown.numerator,mindown.denominator] if mindown is not None else None}
    out=Path(f'/mnt/data/resexact_{k:02d}_{B}_{idx}.json');out.write_text(json.dumps(cert,separators=(',',':')))
    if verify:
        # Recompute all exact conductances from stored original graph weights.
        rows2=[]
        for sr in cert['stresses']:
            nums=list(map(int,sr['weight_nums']));pp,dn,ex=exact_D(B,I,faces,nums)
            assert pp==pairs and dn==list(map(int,sr['Dnums']))
            rows2.append([F(n,QD) for n in dn])
        oldmap={old:(lo,hi) for old,lo,hi in oldleaves};vl=ve=0;vm=None
        for r in roots:
            nleaf,nempty,mm=verify_tree(*oldmap[int(r['old_node'])],r['nodes'],pairs,rows2)
            vl+=nleaf;ve+=nempty;vm=mm if vm is None or (mm is not None and mm<vm) else vm
        assert vl+ve>0
    else:vl=ve=0;vm=None
    log={'entry':k,'B':B,'I':I,'idx':idx,'stresses':len(rows),'old_leaves':len(oldleaves),
         'nodes':tot,'generated_minmargin':None if mn is None else float(mn),
         'verified_leaves':vl,'verified_empty':ve,'verified_minmargin':None if vm is None else float(vm),
         'sec':time.time()-t,'file':str(out)}
    Path(f'/mnt/data/resexact_{k:02d}_{B}_{idx}.log').write_text(json.dumps(log,indent=2))
    print('DONE',json.dumps(log),flush=True)
    return log

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--entry',type=int);ap.add_argument('--start',type=int);ap.add_argument('--end',type=int);ap.add_argument('--tag',default='x');ap.add_argument('--no-verify',action='store_true');a=ap.parse_args()
    inds=[a.entry] if a.entry is not None else list(range(a.start,a.end))
    t=time.time()
    for k in inds:
        try:process_entry(k,not a.no_verify)
        except Exception as ex:
            traceback.print_exc();Path(f'/mnt/data/resexact_{k:02d}.error').write_text(repr(ex))
    print('WORKER_DONE',a.tag,'sec',time.time()-t,flush=True)
