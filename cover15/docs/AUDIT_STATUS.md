# Cover15 audit and comparison

## Input provenance

The two archives and standalone manuscript were received and checked on
September 5, 2026. Their SHA-256 digests are:

| Input | SHA-256 |
|---|---|
| `r15_verified_proof.zip` | `9373b1760afd6890bf4f1c204004ea0796f586b6bcd65a7a10bbc7e0f47e047e` |
| `r15_verified_proof(1).zip` | `9f6bfac36ddc7fd54e2c7d71669885e1e564c9c49921e277dacfc67a0a640eec` |
| `r15_proof_zh.md` | `9a34ccde248e716d4ed128846ca04fb381b39b1643551e3d14bd4f3b246da71c` |

Both archives extracted successfully and their supplied checksum lists passed
before execution. The standalone manuscript is byte-identical to the
alternative archive's `PROOF_zh.md`. The primary manuscript is a different
version and is retained in its own bundle. The repository manifest additionally
covers every archived file, including discovery sources and submitted reports.

## Certificate comparison

| Item | Primary | Alternative |
|---|---:|---:|
| KKT variables | 148 | 148 |
| Geometry coordinates | 76 | 76 |
| Active rods | 61 | 61 |
| Upper complex (V,E,F) | (26,64,39) | (26,64,39) |
| Topology orbits | 11,950,884 | 11,950,884 |
| Metric exclusions | 11,907,870 | 11,907,870 |
| Noncandidate cases | 43,013 | 43,013 |
| Candidate nodes | 983 | 2,029 |
| Candidate force leaves | 459 | 939 |
| Candidate local leaves | 33 | 76 |
| Noncandidate nodes | 58,905 | 58,905 |
| Noncandidate force leaves | 50,959 | 50,959 |
| Anchor radius | 13/100 | 1/9 |

The primary census auditor uses invariant pruning and internal-label
permutations. Its pruning is justified in its manuscript and it additionally
compares selected stabilizers with full group traversal. The alternative uses
breadth-first rooted normalization, distinct from its generator's depth-first
normalization; both use the storage container in `state_store.hpp`. In each
chain, validity, orbit uniqueness and orbit mass are checked against a
separately evaluated rooted counting recurrence. Flip-search termination is
not taken as a completeness proof.

The archived reports are input provenance, not the evidence for a new replay.
Fresh logs and reports are stored separately under `../verification/`.

## Fresh replay results

Both archives were extracted into fresh directories with no cached census.
Both complete master runs returned code 0 and printed
`R15 COMPLETE CERTIFICATE CHAIN VERIFIED`:

| Local replay | Result | Elapsed seconds |
|---|---|---:|
| Primary, `python3 -S -B verify_all.py --jobs 1` | PASS | 575.25 |
| Alternative, `python3 -S -B verify_all.py` | PASS | 108.94 |
| Alternative through the new repository runner | PASS | 204.26 for the inner master |

These timings include regeneration and were obtained under concurrent local
load; they are not a controlled performance comparison. The primary supplied
report was faster but does not replace this fresh regeneration result.

The root, upper cover, local matrix, five census audits, metric partition,
all candidate and noncandidate leaves, and rejection tests passed in each
independent run. Both report zero unresolved leaves. The 146-file repository
manifest passed before and after temporary-copy replay, and `make integrity`
passed for all Cover11 through Cover15 inputs.

The fresh master reports and logs are named
`primary_independent_MASTER_VERIFIED.json`,
`primary_independent_full_replay.log`,
`alternative_independent_MASTER_VERIFIED.json`, and
`alternative_independent_full_replay.log`. The repository entry-point log is
`repository_alternative_replay.log`. Exact root correspondence is recorded in
`bundle_linkage_verified.json`.

Additional tests displaced the alternative radius-squared coordinate and
removed a full face in memory; the new linkage checker rejected both changes
and continued to accept the original inputs. No certificate inconsistency or
specific mathematical error blocking publication was identified in this review.

## Mathematical review

The review covered the strict-radius perturbation to ordinary Voronoi cells,
irreducibility and private points, boundary no-return argument, disk-dual
reduction, stellar padding for fewer essential disks, rooted enumeration
recurrence, safe angle and Farkas constraints, force inequalities, and the
anchor-only isolation theorem. The papers distinguish 148 KKT variables,
76 geometry coordinates and ten independent candidate anchor angles.

The upper proof uses common-point coverage of center triangles, 22 ring
triangles and eleven circular caps. Its disk image covers the anchor polygon
by a boundary winding-number argument. There are thirteen exact active faces
and four strictly slack faces. The primary upper verifier explicitly checks
global vertex connectivity; the repository linkage check explicitly checks
that same property for both upper complexes, supplementing the alternative
verifier's incidence and link checks.

The local certificate is applied to a 61-rod subgraph of the full 43-node
auxiliary graph. The verifiers check the face-set map, subgraph embedding,
all residual indices and complete tree reachability. The proof's constants
and local radii differ between bundles and must not be interchanged.

## Exact root correspondence

The primary-to-alternative center permutation is

```text
[5,6,7,8,9,10,0,1,2,3,4,13,14,11,12].
```

Anchors follow `i -> (i+5) mod 11`; the active witness permutation is checked
from the actual face tables. `verify_bundle_linkage.py` freshly runs both
148-variable root verifiers. It then checks all full faces and active rods,
rotates the alternative anchor 5 to `(1,0)`, relabels coordinates and weights,
and recovers the omitted anchor multiplier using zero total torque. These
operations preserve the KKT equations. Exact interval arithmetic places all
148 transformed coordinates inside the primary `10^-60` uniqueness box.
The root uniqueness theorem therefore identifies the two `t` coordinates;
the conclusion does not rely on matching decimal prefixes alone.

## Limits

The packages share the exact trigonometric and root-verifier source verbatim,
as well as related mathematical arguments. They are related implementations,
not two independent external validations. A successful replay certifies the
finite calculations accepted by these programs; written geometric and
topological reductions remain part of the mathematical review boundary.
Neither peer review nor Lean/Coq formalization is claimed.
