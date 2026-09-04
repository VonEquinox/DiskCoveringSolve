#!/usr/bin/env python3
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
import math,json

# Rigorous pi enclosure by Machin's formula.
def atan_inv_bounds(q:int,nlo:int=31,nhi:int=30):
    def partial(N):
        s=F(0)
        for k in range(N+1):
            z=F(1,(2*k+1)*q**(2*k+1))
            s += z if k%2==0 else -z
        return s
    # odd partial is below, even partial is above
    return partial(nlo),partial(nhi)
a5L,a5U=atan_inv_bounds(5);a239L,a239U=atan_inv_bounds(239)
PI_L=16*a5L-4*a239U
PI_U=16*a5U-4*a239L
S_L=2*PI_L;S_U=2*PI_U
assert PI_L<PI_U
DELTA=F(19489,25000)
TARGET=F(1443877286317993,10**16)

# Fixed-point downward conductance denominator.
QCOEF=1<<90

def floor_frac(x:F,Q:int=QCOEF)->F:
    assert x>=0
    return F((x*Q).numerator//(x*Q).denominator,Q)

from cover11_fixed_trig import Q as TQ, enc as trig_enc, point_nearest, upint, sin_iv, cos_iv

@lru_cache(maxsize=300000)
def line_exact(l:F,u:F):
    assert F(0)<=l<=u<=S_U
    amid=(l+u)/2
    ai=point_nearest(amid); a=F(ai,TQ)
    h=max(a-l,u-a); hi=upint(h); hd=F(hi,TQ)
    sv=sin_iv(__import__('cover11_fixed_trig').I(ai));cv=cos_iv(__import__('cover11_fixed_trig').I(ai))
    si=(sv.lo+sv.hi)//2; s0=F(si,TQ)
    eps=F(max(si-sv.lo,sv.hi-si),TQ)
    fl=F(TQ-cv.hi,TQ)
    if l<=PI_U and u>=PI_L:
        mi=-TQ
    elif u<PI_L:
        mi=cos_iv(trig_enc(u)).lo
    elif l>PI_U:
        mi=cos_iv(trig_enc(l)).lo
    else:
        mi=-TQ
    if mi>0:mi=0
    b=fl-s0*a-eps*hd+F(mi,TQ)*hd*hd/2
    return s0,b

def tighten_exact(lo,hi):
    lo=list(lo);hi=list(hi)
    for _ in range(30):
        slo=sum(lo,F(0));shi=sum(hi,F(0))
        if shi<S_L or slo>S_U:return None
        changed=False
        for i in range(len(lo)):
            nl=S_L-(shi-hi[i]);nh=S_U-(slo-lo[i])
            if nl>lo[i]:lo[i]=nl;changed=True
            if nh<hi[i]:hi[i]=nh;changed=True
            if lo[i]>hi[i]:return None
        if not changed:break
    return tuple(lo),tuple(hi)

def pair_intervals_exact(lo,hi,pairs):
    slo=sum(lo,F(0));shi=sum(hi,F(0));out=[]
    for i,j in pairs:
        l0=sum(lo[i:j],F(0));u0=sum(hi[i:j],F(0))
        l=max(l0,S_L-(shi-u0));u=min(u0,S_U-(slo-l0))
        assert l<=u
        out.append((l,u))
    return out

def linear_min_exact(c,b,lo,hi):
    # minimum over box and S_L <= sum x <= S_U
    n=len(lo);x=[hi[i] if c[i]<0 else lo[i] for i in range(n)]
    v=b+sum((c[i]*x[i] for i in range(n)),F(0));S=sum(x,F(0))
    if S<S_L:
        need=S_L-S
        for i in sorted(range(n),key=lambda z:(c[z],z)):
            d=min(hi[i]-x[i],need);v+=c[i]*d;x[i]+=d;need-=d
            if need==0:break
        assert need==0
    elif S>S_U:
        need=S-S_U
        for i in sorted(range(n),key=lambda z:(c[z],z),reverse=True):
            d=min(x[i]-lo[i],need);v-=c[i]*d;x[i]-=d;need-=d
            if need==0:break
        assert need==0
    return v

def load_D_down(path='/mnt/data/cover11_exact_stresses.json'):
    E=json.load(open(path));pairs=[tuple(map(int,p)) for p in E['pairs']]
    rows=[]
    for row in E['rows']:
        rows.append([floor_frac(F(int(p),int(q))) for p,q in row['D']])
    return pairs,rows

def mixed_row(rows,wit):
    if wit['kind']=='row':return rows[int(wit['row'])]
    nums=list(map(int,wit['nums']));den=int(wit['den']);assert min(nums)>=0 and sum(nums)==den
    out=[]
    for p in range(len(rows[0])):
        z=sum((F(nums[k],den)*rows[k][p] for k in range(len(rows))),F(0))
        out.append(floor_frac(z))
    return out

def bound_exact(lo,hi,wit,pairs,rows):
    th=tighten_exact(lo,hi)
    if th is None:return None
    lo,hi=th;lines=[line_exact(l,u) for l,u in pair_intervals_exact(lo,hi,pairs)]
    row=mixed_row(rows,wit);b=F(0);c=[F(0)]*len(lo)
    for coef,(i,j),(a,b0) in zip(row,pairs,lines):
        b+=coef*b0;da=coef*a
        for h in range(i,j):c[h]+=da
    return linear_min_exact(c,b,lo,hi)
