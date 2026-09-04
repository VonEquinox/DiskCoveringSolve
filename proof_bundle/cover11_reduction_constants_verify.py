#!/usr/bin/env python3
"""Exact rational/Taylor checks for the global reduction constants.

These checks certify that every boundary gap in a hypothetical cover with
r^2<T+ lies below DELTA=0.77956, that eight such arcs cannot cover the unit
circle, and that r<1/2.
"""
from fractions import Fraction as F
import sys,json
sys.path.insert(0,'/mnt/data')
import cover11_candidate_exact_core as core
from cover11_fixed_trig import Q,enc,sin_iv
T=core.TARGET
DELTA=core.DELTA
s=sin_iv(enc(DELTA/2)); slo=F(s.lo,Q)
# If an arc has angle alpha>=DELTA and alpha<pi, its endpoint chord is at
# least 2 sin(DELTA/2), hence cannot lie in one radius-r disk when r^2<T.
assert slo>0
arc_margin=slo*slo-T
assert arc_margin>0
# At most eight boundary cells is impossible.
length_margin=core.S_L-8*DELTA
assert length_margin>0
# A disk whose Voronoi cell contains the origin cannot also touch the unit
# boundary: 2r<1.
half_margin=1-4*T
assert half_margin>0
out={
 'T':[T.numerator,T.denominator],
 'DELTA':[DELTA.numerator,DELTA.denominator],
 'sin_delta_over_2_lower':[slo.numerator,slo.denominator],
 'sin2_minus_T':[arc_margin.numerator,arc_margin.denominator],
 'two_pi_lower_minus_8delta':[length_margin.numerator,length_margin.denominator],
 'one_minus_4T':[half_margin.numerator,half_margin.denominator],
 'floats':{
  'sin2_minus_T':float(arc_margin),
  'two_pi_lower_minus_8delta':float(length_margin),
  'one_minus_4T':float(half_margin),
 }
}
open('/mnt/data/cover11_reduction_constants.json','w').write(json.dumps(out,indent=2))
print('GLOBAL REDUCTION CONSTANTS VERIFIED')
for k,v in out['floats'].items():print(k,v)
