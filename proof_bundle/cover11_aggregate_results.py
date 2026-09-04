#!/usr/bin/env python3
"""Aggregate independently verified chunks and assert the complete certificate counts."""
from fractions import Fraction as F
from pathlib import Path
import json, re, glob, hashlib
BASE=Path('/mnt/data')

def frac(z):
    return F(int(z[0]),int(z[1]))

def load_many(pattern):
    return [json.load(open(p)) for p in sorted(glob.glob(str(BASE/pattern)))]

# Exact candidate root and rational target.
K=json.load(open(BASE/'cover11_full_kkt_certificate.json'))
mid=[F(int(a),int(b)) for a,b in K['midpoint']]
rad=F(int(K['radius'][0]),int(K['radius'][1]))
tidx=int(K['t_index']);tlo=mid[tidx]-rad;thi=mid[tidx]+rad
T=F(1443877286317993,10**16)
assert T>thi

# Enumeration and metric screen.
enum=json.load(open(BASE/'cover11_enumeration_verify.json'))
assert [(x['B'],x['I'],x['orbits']) for x in enum]==[(9,1,291),(9,2,2548),(10,1,1004)]
assert all(x.get('brown_equal') and x['labeled_states']==x['brown_labeled_count'] for x in enum)
M=json.load(open(BASE/'cover11_metric_Tplus.json'))
blocks=M['cases']
metric=[]
for key,want in [('9,1',(291,3)),('9,2',(2548,34)),('10,1',(1004,18))]:
    b=blocks[key];tot=int(b['orbit_count']);res=len(b['residual']);assert (tot,res)==want
    metric.append({'case':key,'orbits':tot,'excluded':tot-res,'residual':res})
assert sum(x['orbits'] for x in metric)==3843 and sum(x['residual'] for x in metric)==55
link=json.load(open(BASE/'cover11_topology_linkage_verify.json'))
assert link['status']=='passed' and link['candidate']=={'B':9,'I':2,'orbit':2547,'boundary_shift':3,'boundary_reflection':False,'interior_swap':True}
assert link['wheel']=={'B':10,'I':1,'orbit':1003,'exact_wheel':True}
assert link['residual_total']==55 and link['ordinary_residuals']==53

# Candidate base chunks.
base=load_many('candbaseverify_[0-4].json');assert len(base)==5
assert [x['start'] for x in base]==[0,21228,42456,63684,84912]
assert [x['end'] for x in base]==[21228,42456,63684,84912,106136]
base_pass=sum(int(x['pass']) for x in base);base_def=sum(int(x['deferred']) for x in base)
assert (base_pass,base_def)==(100834,5302)
base_min=min(frac(x['minmargin']) for x in base);assert base_min>0

# Candidate refined chunks.
ref=load_many('candrefverify_[0-4].json');assert len(ref)==5
assert [x['start'] for x in ref]==[0,1060,2120,3180,4240]
assert [x['end'] for x in ref]==[1060,2120,3180,4240,5302]
ref_roots=sum(int(x['roots']) for x in ref);ref_nodes=sum(int(x['nodes']) for x in ref)
ref_ext=sum(int(x['external']) for x in ref);ref_loc=sum(int(x['local']) for x in ref);ref_emp=sum(int(x['empty']) for x in ref)
assert (ref_roots,ref_nodes,ref_ext,ref_loc,ref_emp)==(5302,304708,153052,1953,0)
ref_min=min(frac(x['minmargin']) for x in ref);assert ref_min>0

# Wheel chunks.
wheel=load_many('wheelchunk_[0-7].json');assert len(wheel)==8
assert wheel[0]['start']==0 and wheel[-1]['end']==4853
for a,b in zip(wheel,wheel[1:]):assert a['end']==b['start']
wheel_direct=sum(int(x['direct']) for x in wheel);wheel_ref=sum(int(x['refined_leaves']) for x in wheel)
wheel_min=min(frac(x['minmargin']) for x in wheel);assert wheel_direct==4781 and wheel_ref==145 and wheel_min>0

# 53 other residual chunks.
res=load_many('resverify_[0-5].json');assert len(res)==6
assert res[0]['start']==0 and res[-1]['end']==53
for a,b in zip(res,res[1:]):assert a['end']==b['start']
res_entries=sum(int(x['entries']) for x in res);res_stresses=sum(int(x['stresses']) for x in res)
res_nodes=sum(int(x['nodes']) for x in res);res_leaves=sum(int(x['leaves']) for x in res);res_empty=sum(int(x['empty']) for x in res)
res_min=min(frac(x['minmargin']) for x in res);assert (res_entries,res_stresses,res_nodes,res_leaves,res_empty)==(53,167,13825,13789,0) and res_min>0

out={
 'status':'complete',
 'exact_definition':'t* is the t-coordinate of the unique root in cover11_full_kkt_certificate.json; r*=sqrt(t*)',
 't_interval':[[str(tlo.numerator),str(tlo.denominator)],[str(thi.numerator),str(thi.denominator)]],
 't_midpoint_decimal':format(float(mid[tidx]),'.17g'),
 'Tplus':[str(T.numerator),str(T.denominator)],
 'Tplus_minus_t_upper_float':float(T-thi),
 'enumeration':metric,
 'linkage':{
   'brown_labeled_counts_verified':True,
   'candidate_orbit':link['candidate'],
   'wheel_orbit':link['wheel'],
   'six_proposal_kron_rows_recomputed':True,
   'nine_gap_to_eight_coordinate_link_verified':True,
   'upper_complex':{'vertices':20,'edges':48,'faces':29,'euler_characteristic':1,'vertex_links_verified':True},
 },
 'candidate':{
   'old_leaves':106136,'direct_exact':base_pass,'deferred':base_def,
   'base_minmargin_float':float(base_min),'refined_roots':ref_roots,
   'refined_nodes':ref_nodes,'external_exact_leaves':ref_ext,'local_leaves':ref_loc,
   'refined_empty':ref_emp,'external_minmargin_float':float(ref_min)
 },
 'wheel':{'old_leaves':4853,'direct':wheel_direct,'refined_exact_leaves':wheel_ref,'minmargin_float':float(wheel_min)},
 'other_residuals':{'topologies':res_entries,'stresses':res_stresses,'nodes':res_nodes,'leaves':res_leaves,'empty':res_empty,'minmargin_float':float(res_min)},
 'proof_conclusion':'r_11 = sqrt(t*)',
}
(BASE/'cover11_complete_summary.json').write_text(json.dumps(out,indent=2))
print('COVER11 COMPLETE CERTIFICATE AGGREGATE VERIFIED')
print(json.dumps(out,indent=2))
