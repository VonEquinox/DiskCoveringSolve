#!/usr/bin/env python3
"""Exact structural audit for the candidate topology certificate partition.

Checks that:
* the 13 stored stress rows come from nonnegative unit-sum graph weights and
  exact star-mesh elimination;
* the five base-witness chunks cover every old leaf exactly once;
* every base-deferred leaf has exactly one refined root;
* refined trees are genuine binary partitions of their old boxes.
Numerical inequalities at leaves are replayed by the separate base/refined
Fraction verifiers.
"""
from fractions import Fraction as F
import json,glob,sys
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core

E=json.load(open('/mnt/data/cover11_exact_stresses.json'))
edges=[tuple(map(int,e)) for e in E['edges']]
pairs=[tuple(map(int,p)) for p in E['pairs']]
assert pairs==[(i,j) for i in range(9) for j in range(i+1,9)]
assert len(edges)==45 and len(E['rows'])==13

def exact_kron(nums,den):
    adj={i:{} for i in range(29)}
    def add(a,b,z):
        if a==b or z==0:return
        adj[a][b]=adj[a].get(b,F(0))+z
        adj[b][a]=adj[b].get(a,F(0))+z
    for (a,b),n in zip(edges,nums):add(a,b,F(n,den))
    for v in range(9,29):
        nei=[(a,z) for a,z in adj.get(v,{}).items() if z]
        s=sum((z for _,z in nei),F(0))
        if s:
            for k,(a,x) in enumerate(nei):
                for b,y in nei[k+1:]:add(a,b,x*y/s)
            for a,_ in nei:adj[a].pop(v,None)
        adj.pop(v,None)
    return [2*adj.get(i,{}).get(j,F(0)) for i,j in pairs]

for k,row in enumerate(E['rows']):
    nums=list(map(int,row['weight_nums']));den=int(row['weight_den'])
    assert len(nums)==45 and min(nums)>=0 and sum(nums)==den
    D=[F(int(a),int(b)) for a,b in row['D']]
    assert D==exact_kron(nums,den),(k,'Kron mismatch')

old=json.load(open('/mnt/data/cover11_candidate_cert_T14438.json'))['nodes']
stack=[(0,(F(0),)*9,(F(1),)*9)];oldboxes={};leaf_order=[];seen=set()
while stack:
    idx,lo,hi=stack.pop();assert idx not in seen;seen.add(idx);n=old[idx]
    if 'split' in n:
        k=int(n['split']);m=(lo[k]+hi[k])/2
        h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
        stack.append((int(n['child'][1]),tuple(l1),hi))
        stack.append((int(n['child'][0]),lo,tuple(h0)))
    elif 'leaf' in n:
        oldboxes[idx]=(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi));leaf_order.append(idx)
    elif 'empty' in n:
        assert core.tighten_exact(tuple(x*core.DELTA for x in lo),tuple(x*core.DELTA for x in hi)) is None
    else:raise AssertionError((idx,n))
assert len(seen)==len(old) and len(leaf_order)==106136

base=[]
for f in sorted(glob.glob('/mnt/data/basewit_[0-4].json')):base+=json.load(open(f))['records']
base.sort(key=lambda r:int(r['rank']))
assert [int(r['rank']) for r in base]==list(range(len(leaf_order)))
assert [int(r['node']) for r in base]==leaf_order
assert len(base)==106136
deferred={int(r['node']) for r in base if r['witness'] is None or r['witness'].get('kind')=='fail'}
direct={int(r['node']) for r in base}-deferred
assert len(deferred)==5302 and len(direct)==100834

R=json.load(open('/mnt/data/cover11_candidate_refined_rational_splits.json'))['roots']
rmap={int(r['old_node']):r for r in R}
assert len(rmap)==len(R)==5303
assert deferred<=set(rmap)
extra=set(rmap)-deferred
assert extra=={20052} and 20052 in direct

counts={'external':0,'local':0,'empty':0,'split':0}
for oldnode in sorted(deferred):
    r=rmap[oldnode];nodes=r['nodes'];boxes=[None]*len(nodes);boxes[0]=oldboxes[oldnode]
    for idx,n in enumerate(nodes):
        assert boxes[idx] is not None,(oldnode,idx)
        lo,hi=boxes[idx]
        if 'split' in n:
            k=int(n['split']);m=F(int(n['mid'][0]),int(n['mid'][1]));assert lo[k]<m<hi[k]
            h0=list(hi);h0[k]=m;l1=list(lo);l1[k]=m
            c0,c1=map(int,n['child']);assert boxes[c0] is None and boxes[c1] is None
            boxes[c0]=(lo,tuple(h0));boxes[c1]=(tuple(l1),hi);counts['split']+=1
        elif n.get('kind') in ('row','mix'):counts['external']+=1
        elif 'local' in n:counts['local']+=1
        elif 'empty' in n:counts['empty']+=1
        else:raise AssertionError((oldnode,idx,n))
    assert all(b is not None for b in boxes)
assert counts['external']==153052 and counts['local']==1953 and counts['empty']==0
print('CANDIDATE CERTIFICATE STRUCTURE VERIFIED')
print('direct',len(direct),'deferred',len(deferred),'refined_extra',sorted(extra))
print(counts)
