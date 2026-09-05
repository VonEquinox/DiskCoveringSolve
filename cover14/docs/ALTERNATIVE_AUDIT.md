# Alternative Cover14 package comparison

The archive `r14_verified_proof(1).zip` supplies a useful alternative certificate
chain for the same algebraic candidate. Its archive SHA-256 is
`6dd4f1f40bdc76a04e32eaae47a59b898bd772afe12f37838d3e9f49ff14f09c`.
All 64 entries in its original manifest passed before replay. Its exact full
replay on September 5, 2026 returned successfully in approximately 35 seconds.
The fresh master and transcript are archived in `../verification/` as
`alternative_independent_MASTER_VERIFIED.json` and
`alternative_independent_full_replay.log`.

## Comparison

| Item | Primary package | Alternative package |
|---|---:|---:|
| KKT variables | 116 | 116 |
| Topology orbits | 1,313,024 | 1,313,024 |
| Metric exclusions | 1,305,792 | 1,305,792 |
| Noncandidate cases | 7,231 | 7,231 |
| Candidate tree nodes | 913 | 599 |
| Candidate force leaves | 419 | 266 |
| Candidate local leaves | 38 | 34 |
| Noncandidate tree nodes | 16,075 | 11,621 |
| Noncandidate force leaves | 11,653 | 9,426 |
| Anchor isolation radius | 1/6 | 7/40 |
| Upper complex (V,E,F) | (24,59,36) | (24,57,34) |

Both replays have no unresolved cases and certify exactly the same displayed
strict intervals for the radius and its reciprocal. The primary upper proof
allows two collapsed triangles in a continuous disk map. The alternative
removes those triangles and inserts their centers into the boundary chords,
leaving 34 strictly positively oriented triangles. Both written arguments
account for the degeneracy.

The alternative manuscript also explicitly describes the one-index anchor
offset under a reflection of boundary-center labels. Its discovery directory
includes additional numeric checkpoints and generation scripts.

## Exact linkage

`../verify_bundle_linkage.py` checks the full face sets, active witness sets,
and all 44 rods under the primary-to-alternative center permutation

```text
[3,4,5,6,7,8,9,0,1,2,11,12,13,10].
```

The corresponding anchor map is `i -> (i+3) mod 10`. The alternative root is
rotated so its anchor 3 becomes `(1,0)`, all points and weights are relabeled,
and the omitted anchor multiplier is recovered from the radial force balance.
The total-torque identity justifies that recovery, as in the symmetry proofs.

The graph map preserves squared rod lengths, unit-anchor equations and KKT
stationarity under this isometry. Exact rational interval arithmetic places
the resulting 116 coordinates strictly inside the primary uniqueness box of
radius `10^-60`. Together with the root existence and uniqueness checks in the
two full replays, this identifies the algebraic roots, including their `t`
coordinate. Matching decimal prefixes alone is not used for this conclusion.

## Reproduction and scope

```bash
make cover14-alternative-full PYTHON_BIN="$(command -v python3)"
```

The repository runner checks the original manifest, copies the archived bundle
to a temporary directory, runs the entire verifier, checks the root linkage,
and rechecks the archived manifest. Original submitted reports remain fixed.
The primary replay remains available through `make cover14-full`.

These are related certificate implementations, not evidence of two independent
research teams or external peer review. In particular, `exact_arcs.py` is
byte-identical between the two archives, and their orbit audit logs agree.
Their enumeration files differ even after decompression, so residual indices
and tree data must stay with their own verifier. Each package separately
recomputes validity, orbit uniqueness and completeness. The alternative has
smaller trees and a larger certified anchor domain, but retains shared
mathematical and software dependencies.
