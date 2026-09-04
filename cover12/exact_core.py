#!/usr/bin/env python3
"""Small exact-arithmetic kernel for congruent-disk covering certificates.

All transcendental inequalities are reduced to rational Taylor bounds.  No
floating-point result is trusted by the verification routines.
"""
from fractions import Fraction as F
from math import factorial


def atan_inv_bounds(q: int, terms: int = 70):
    """Rational enclosure of atan(1/q) from its alternating series."""
    assert q >= 2 and terms >= 1
    s = F(0)
    x = F(1, q)
    x2 = x*x
    p = x
    for k in range(terms):
        term = p / (2*k+1)
        s += term if (k & 1) == 0 else -term
        p *= x2
    nxt = p / (2*terms+1)
    # terms counts k=0,...,terms-1; next sign is (-1)^terms.
    if terms & 1:
        # next term negative: partial sum is an upper bound.
        return s-nxt, s
    else:
        return s, s+nxt


def machin_pi_bounds(terms5: int = 70, terms239: int = 24):
    a5l,a5u = atan_inv_bounds(5,terms5)
    a239l,a239u = atan_inv_bounds(239,terms239)
    # pi = 16 atan(1/5) - 4 atan(1/239)
    return 16*a5l-4*a239u, 16*a5u-4*a239l

PI_L, PI_U = machin_pi_bounds()
TWO_PI_L, TWO_PI_U = 2*PI_L, 2*PI_U


def sincos_interval(x: F, n: int = 42):
    """Return rational intervals containing sin(x), cos(x).

    Taylor's theorem with |derivative|<=1 is used, so this is valid for every
    rational x.  The intended range is |x|<=2*pi.
    """
    x = F(x)
    ax = abs(x)
    # sin through degree 2n+1; first omitted Taylor degree is 2n+2, but that
    # coefficient is zero.  We use the general degree-(2n+1) remainder bound.
    ss = F(0)
    term = x
    for k in range(n+1):
        if k:
            term *= -x*x / ((2*k)*(2*k+1))
        ss += term
    rs = ax**(2*n+2) / factorial(2*n+2)

    cc = F(0)
    term = F(1)
    for k in range(n+1):
        if k:
            term *= -x*x / ((2*k-1)*(2*k))
        cc += term
    rc = ax**(2*n+1) / factorial(2*n+1)
    return ss-rs, ss+rs, cc-rc, cc+rc


def phi_interval(x: F, n: int = 42):
    sl,su,cl,cu = sincos_interval(x,n)
    return 1-cu,1-cl


def cos_lower_on_interval(l: F, u: F, n: int = 42):
    """A rational lower bound for min_{[l,u]} cos on 0<=l<=u<=2*pi.

    Uses monotonicity on [0,pi] and [pi,2pi], with rigorous pi bounds.  If the
    interval might straddle pi, the universal lower bound -1 is returned.
    """
    assert 0 <= l <= u <= TWO_PI_U
    if u <= PI_L:
        return sincos_interval(u,n)[2]
    if l >= PI_U:
        return sincos_interval(l,n)[2]
    return F(-1)


def line_exact(l: F, u: F, n: int = 42):
    """Rational affine minorant A*x+B <= 1-cos(x) on [l,u]."""
    l,u=F(l),F(u)
    assert l <= u
    a=(l+u)/2
    h=(u-l)/2
    sl,su,cl,cu=sincos_interval(a,n)
    A=(sl+su)/2
    eps=(su-sl)/2
    phil=1-cu
    mc=cos_lower_on_interval(l,u,n)
    # Taylor with phi''=cos.  If mc>0 the quadratic term is nonnegative;
    # otherwise its minimum over |x-a|<=h is mc*h^2/2.
    q = min(F(0),mc)*h*h/2
    B=phil-A*a-eps*h+q
    return A,B


def tighten_sum_box(lo,hi,L=TWO_PI_L,U=TWO_PI_U):
    lo=list(map(F,lo));hi=list(map(F,hi))
    for _ in range(64):
        slo=sum(lo,F(0));shi=sum(hi,F(0))
        if shi < L or slo > U:return None
        changed=False
        for i in range(len(lo)):
            nl=L-(shi-hi[i]);nh=U-(slo-lo[i])
            if nl>lo[i]:lo[i]=nl;changed=True
            if nh<hi[i]:hi[i]=nh;changed=True
            if lo[i]>hi[i]:return None
        if not changed:break
    return tuple(lo),tuple(hi)


def pair_intervals(lo,hi,pairs,L=TWO_PI_L,U=TWO_PI_U):
    slo=sum(lo,F(0));shi=sum(hi,F(0));out=[]
    for i,j in pairs:
        ml0=sum(lo[i:j],F(0));mh0=sum(hi[i:j],F(0))
        ml=max(ml0,L-(shi-mh0));mh=min(mh0,U-(slo-ml0))
        assert ml<=mh
        out.append((ml,mh))
    return out


def linear_min(c,lo,hi,L=TWO_PI_L,U=TWO_PI_U):
    """Exact min c.x on box intersected with L<=sum x<=U."""
    n=len(c);x=[hi[i] if c[i]<0 else lo[i] for i in range(n)]
    v=sum((c[i]*x[i] for i in range(n)),F(0));S=sum(x,F(0))
    if S<L:
        need=L-S
        for i in sorted(range(n),key=lambda z:c[z]):
            d=min(hi[i]-x[i],need);v+=c[i]*d;x[i]+=d;need-=d
            if need==0:break
        assert need==0
    elif S>U:
        need=S-U
        for i in sorted(range(n),key=lambda z:c[z],reverse=True):
            d=min(x[i]-lo[i],need);v-=c[i]*d;x[i]-=d;need-=d
            if need==0:break
        assert need==0
    return v


def graph_edges(B,I,faces,boundary_order=None):
    """Labels and edges of the anchor/center/face witness graph."""
    if boundary_order is None:boundary_order=list(range(B))
    anchors=[f'q{i}' for i in range(B)]
    centers=[f'c{i}' for i in range(B+I)]
    face_nodes=[f'p{k}' for k in range(len(faces))]
    edges=[]
    for i in range(B):
        edges += [(f'q{i}',f'c{boundary_order[i]}'),
                  (f'q{i}',f'c{boundary_order[(i+1)%B]}')]
    for k,fa in enumerate(faces):
        edges += [(f'p{k}',f'c{v}') for v in fa]
    return anchors,centers+face_nodes,edges


def exact_kron(B,I,faces,nums,den,boundary_order=None):
    """Component-safe exact star-mesh elimination.

    Returns coefficients D_ij in E=sum D_ij(1-cos(theta_j-theta_i)).
    The input edge weights are nums/den and must be nonnegative, summing to 1.
    Free-only components are irrelevant and discarded.
    """
    nums=list(map(int,nums));den=int(den)
    A,free,E=graph_edges(B,I,faces,boundary_order)
    assert len(nums)==len(E) and min(nums)>=0 and sum(nums)==den
    allv=A+free
    adj={v:{} for v in allv}
    def add(a,b,w):
        if a==b or w==0:return
        adj[a][b]=adj[a].get(b,F(0))+w
        adj[b][a]=adj[b].get(a,F(0))+w
    for (a,b),z in zip(E,nums):add(a,b,F(z,den))

    # Discard components containing no anchor.
    seen=set();keep=set(A)
    for s in allv:
        if s in seen:continue
        stack=[s];seen.add(s);comp=[];has_anchor=False
        while stack:
            v=stack.pop();comp.append(v);has_anchor |= v in A
            for w in adj[v]:
                if w not in seen:seen.add(w);stack.append(w)
        if has_anchor:keep.update(comp)
    for v in list(adj):
        if v not in keep:
            for w in list(adj[v]):adj[w].pop(v,None)
            del adj[v]

    for v in free:
        if v not in adj:continue
        nei=list(adj[v].items());S=sum((w for _,w in nei),F(0))
        if S==0:
            del adj[v];continue
        for x in range(len(nei)):
            a,wa=nei[x]
            for y in range(x+1,len(nei)):
                b,wb=nei[y];add(a,b,wa*wb/S)
        for a,_ in nei:adj[a].pop(v,None)
        del adj[v]
    pairs=[(i,j) for i in range(B) for j in range(i+1,B)]
    # Conductance c_ij multiplies |q_i-q_j|^2=2(1-cos), hence D=2c.
    return [2*adj[A[i]].get(A[j],F(0)) for i,j in pairs]


if __name__=='__main__':
    print('pi width',float(PI_U-PI_L))
    for z in [F(0),F(1,10),F(1),F(3),F(6)]:
        sl,su,cl,cu=sincos_interval(z)
        assert sl<=su and cl<=cu
    print('self-test ok')

# --- Tighter affine minorants used by the final certificate ---
from functools import lru_cache
from math import isqrt, asin

@lru_cache(maxsize=200000)
def sincos_cached(x, n=34):
    return sincos_interval(F(x),n)

def sqrt_upper(y: F, qbits: int = 100):
    """Dyadic rational upper bound for sqrt(y), y>=0."""
    y=F(y);assert y>=0
    Q=1<<qbits
    # N=ceil(y Q^2), then ceil(sqrt(N))/Q is certainly >=sqrt(y).
    z=y*Q*Q;N=z.numerator//z.denominator + (z.numerator%z.denominator!=0)
    n=isqrt(N)
    if n*n<N:n+=1
    return F(n,Q)

@lru_cache(maxsize=100000)
def asin_bracket_nonnegative(A: F, qbits: int = 46):
    """Verified rational bracket for asin(A), 0<=A<1."""
    A=F(A);assert 0<=A<1
    Q=1<<qbits
    z=asin(float(A));n=int(z*Q)
    # Start with a few ulps around the floating proposal and enlarge until the
    # exact Taylor signs verify the bracket.
    pad=2
    while True:
        lo=F(max(0,n-pad),Q);hi=F(n+pad+1,Q)
        # Keep the upper endpoint in the monotone interval [0,pi/2].
        if 2*hi>=PI_L:
            hi=PI_L/2
        slo=sincos_cached(lo,34)[0]
        shi=sincos_cached(hi,34)[1]
        if slo<=A<=shi:return lo,hi
        pad*=2
        assert pad<1<<20

@lru_cache(maxsize=200000)
def line_exact_tight(l: F,u: F,grid_bits: int = 12):
    """Tight rational affine minorant A*x+B <= 1-cos x on [l,u].

    A is a rational enclosure midpoint of sin(a) at a nearby dyadic grid
    point.  The minimum of 1-cos(x)-A*x occurs at an endpoint or at its one
    local minimum; the latter is enclosed via a verified asin bracket.
    """
    l,u=F(l),F(u);assert 0<=l<=u<=TWO_PI_U
    mid=(l+u)/2;Q=1<<grid_bits
    z=mid*Q;j=(2*z.numerator+z.denominator)//(2*z.denominator)
    a=F(j,Q)
    sl,su,cl,cu=sincos_cached(a,34);A=(sl+su)/2
    # Numerical noise cannot put this rational midpoint outside [-1,1], but
    # assert it exactly.
    assert -1 < A < 1
    candidates=[]
    for x in (l,u):
        phil=1-sincos_cached(x,34)[3]
        candidates.append(phil-A*x)
    if A>0:
        xl,xu=asin_bracket_nonnegative(A)
        if not (xu<l or xl>u):
            root_lo=max(xl,l);root_hi=min(xu,u)
            # If clipping reaches an endpoint, endpoint candidates already
            # handle it.  The following remains a valid extra lower candidate.
            suroot=sqrt_upper(1-A*A)
            candidates.append(1-suroot-A*root_hi)
    elif A<0:
        bl,bu=asin_bracket_nonnegative(-A)
        xl=TWO_PI_L-bu;xu=TWO_PI_U-bl
        if not (xu<l or xl>u):
            root_lo=max(xl,l);root_hi=min(xu,u)
            suroot=sqrt_upper(1-A*A)
            candidates.append(1-suroot-A*root_lo)
    return A,min(candidates)
