# Cover18 audit and release status

Reviewed and replayed locally on September 5, 2026. Both submitted bundles
passed actual full execution, then each passed a second full execution through
the repository runner. No blocking mathematical or exact-arithmetic defect was
identified in this review. This is not a claim of external peer review,
proof-assistant formalization, or exhaustive testing of every software defect.

## Finding and mitigation

The alternative upstream `verify_all.py:11` rejects `-O` before removing its
previous `reports/MASTER_VERIFIED.json`. A failed optimized invocation can
therefore leave an old success report, despite the broad cleanup wording in
the upstream documentation. The process correctly returns nonzero; this does
not affect the successful normal replays or change any certificate inequality.
The primary entry point deletes its old master before checking this option.

This behavior was reproduced with a synthetic old report in a private temporary
directory. Evidence is in `verification/upstream_guard_regression.json`.
The repository wrapper deletes all reports and generated core results in a
disposable copy, strips optimization/PYTHONPATH environment overrides, checks
the subprocess exit code and every expected stage, and exports only new results
after success. It does not rewrite either original entry point or certificate.
Use the repository wrapper rather than treating any supplied JSON flag as a
fresh verification result.

## Original inputs

The primary ZIP contains 31 files, 34,009,699 uncompressed bytes. The supplied
standalone `r18_proof_zh.md` matches its `PROOF_zh.md` byte for byte.
The alternative ZIP contains 76 files, 88,862,765 uncompressed bytes. Its own
75-entry `SHA256SUMS` passed before execution. ZIP CRCs, paths, duplicate names
and symlink checks passed for both archives. The repository manifest preserves
all 107 original files, including the alternative's submitted reports.

| Input | SHA-256 |
| --- | --- |
| Primary ZIP | `e35283e2d9fe3f87d5340269de94dd832c33b87f3ed7c6b0e99afbaa918de729` |
| Alternative ZIP | `cc212ac381af41273ecfff94baf22e65e3de80afdac87384516f86fd567d5957` |
| Standalone primary proof | `b57e1c7c475d6756fc452f4610675320701e4f5eae4dfc926cc1f08e834b8789` |
| Primary root | `446d47409d67bd6d6d2af39813b6c4b974f29d896ea8b25f7ec15bf78dd13c34` |
| Alternative root | `b77a77336488a0f85981df24ddfb897a048113e4cce55453f0175b99c1274911` |

## Actual full replays

Python 3.14.6 at `/opt/homebrew/bin/python3`, standard-library-only execution
with `-S -B`, and the local C++17 compiler were used. Full environment details
are in `verification/environment_and_inputs.json`.

| Run | Primary seconds | Alternative seconds | Result |
| --- | ---: | ---: | --- |
| Fresh independent extraction | 53.9806 | 69.4624 | Both exit 0 |
| Repository runner, disposable copies | 49.0411 | 61.3766 | Both exit 0 |

Times are the individual master-replay durations and exclude the repository's
additional cross-bundle comparison. Stored upstream reports were not accepted
in place of these runs. Original acceptance input hashes remained unchanged.
The second pair additionally checked the repository wrapper, full stage lists,
report export and source preservation end to end.

Each primary run executed all seven accepting stages. Each alternative run
executed all 24 stages, including fresh compiler invocations and separate DFS
and BFS runs for all seven boundary families. Both chains returned zero
unresolved leaves. Their local matrix margins were respectively `793/56250000`
and `6977/78125000`, strictly positive rational numbers.

| Boundary / interior vertices | Rooted survivors | Orbit survivors |
| --- | ---: | ---: |
| 11 / 7 | 8,470 | 397 |
| 12 / 6 | 165,816 | 6,971 |
| 13 / 5 | 186,121 | 7,191 |
| 14 / 4 | 44,926 | 1,645 |
| 15 / 3 | 320 | 16 |
| 16 / 2 | 0 | 0 |
| 17 / 1 | 0 | 0 |
| Total | 405,653 | 16,220 |

The primary visited 3,180,122 partial states and its separate auditor visited
10,018,169. These are not counts of all unpruned topologies. Completeness
depends on the written root-face recursion and safe necessary-condition cuts.

The primary passed 26 unpruned count cases, including `B=3, I=4,5,6,7`; the
alternative passed 29 cases. Each passed 32 deliberate rejection tests. The
alternative also checked three exact trigonometric identities. Twelve new
synthetic wrapper tests passed, including corrupt/missing/extra files, unsafe
manifest paths, mismatched proofs, stale reports, skipped stages, unresolved
results, missing regression checks, export/cleanup and protected paths.
`make integrity` passed for Covers 11-18; this does not claim fresh full replays
of all older covers in this release.

## Exact cross-bundle linkage

`verify_bundle_linkage.py` searches boundary dihedral maps and permutations of
the six interior centers, requiring both complete 22-face tables and active
16-face tables to match. One orientation-preserving primary-to-alternative
center map is:

```text
[2,3,4,5,6,7,8,9,10,11,0,1,17,16,15,14,13,12]
```

The anchors shift by two, and all 72 rods and 16 active witnesses are explicitly
matched. Both root certificates are re-executed before transportation. The
alternative `10^-80` root box, after rotation, node and multiplier permutation,
maps strictly inside the primary `10^-60` uniqueness box. The maximum interval
deviation is less than `2.771e-80`; its exact fraction is recorded in the JSON.

The omitted multiplier at the fixed anchor is recovered from its incident
forces. Stationarity at all other nodes and total torque cancellation make
that remaining force radial. Orthogonal transport therefore preserves the
full constrained stationarity equations, including the newly fixed anchor's
normalization, and leaves `t` unchanged. Uniqueness in the primary box proves
that both isolated roots have the same exact `t`, not merely matching decimals.
A reflection-based transport passes as well but is not needed for equality.

All 16,220 survivor representatives are separately converted to the same
Python BFS canonical form, then their complete sets and stabilizers compared
family by family. This comparison passes despite different internal labels,
candidate row indices, pruning constants and branch-tree data. The packages
must still remain separate because their individual certificate interfaces
are not interchangeable.

## Claim-evidence review

| Claim | Evidence reviewed and replayed | Status |
| --- | --- | --- |
| A unique specified algebraic root exists | Explicit quadratic system; rational contraction on nested boxes | Supported in the specified boxes, not globally unique |
| Eighteen disks cover the entire unit disk | 16 active and 6 strict center-face witnesses; 46 positive triangles; boundary chain/winding argument; 12 circular caps | Supported |
| Every smaller cover enters the finite reduction | Buffered ordinary Voronoi construction, private-point no-return argument, simple dual, boundary bound, padding | Written mathematical argument reviewed; not formalized |
| Partial pruning loses no possible counterexample | Root-face induction and persistent-edge path inequalities; two fresh enumerators per bundle | Supported by written proof plus finite execution |
| All surviving continuous cases are excluded or locally isolated | Exact candidate graph maps, force balances, full tree coverage and positive-definite matrix certificates | Supported; zero unresolved leaves |
| Two archives give the same exact answer | Fresh root-box transport and full orbit-set comparison | Supported |
| Numerical discovery was independently reproduced | Optional proposal sources are included only in the alternative | Not claimed; discovery was not rerun |
| External peer review or Lean/Coq certification is complete | No such evidence supplied | Not claimed |

Both full Chinese proofs were read, including their reference sections. Their
logical order was checked from definition and isolated root through upper
cover, geometric reduction, exhaustive recursion, force/local lemmas and final
contradiction. In particular, neither assumes Covers 11-17, twelve boundary
cells for competitors, candidate symmetry, or local KKT optimality as a global
lower bound. The upper-cover winding argument makes a separate segment-crossing
guess unnecessary. General topology and analysis arguments remain human-readable
proof obligations, not facts established merely by a `PASS` log.

For historical attribution, Friedman's [Circles Covering Circles](https://erich-friedman.github.io/packing/circovcir/)
was checked on the review date: its `n=18`, reciprocal-radius `3.446+` entry
attributes the construction to Jeremy Tan (2018). This agrees with both
proofs' background statements and is not used to justify the lower bound or
any priority claim for this repository.

The review dimensions are claim support, exposition, computational coverage,
evaluation integrity and method justification. Machine-learning baselines,
ablations and statistical significance are not applicable to exact certificate
acceptance. No additional idea report, experiment dev log or evaluation-results
directory was supplied for this task.

## Evidence and maintenance

`verification/independent_primary/` and `independent_alternative/` contain fresh
masters, stage logs and generated results. The corresponding `.log` files hold
actual standard output. `verification/repository_replay/` holds the second pair
and its fresh linkage result. `independent_linkage_verified.json` records the
first extra linkage check. `repository_full_replay.log` and
`integrity_and_runner_tests.log` preserve the integration checks.

CI configuration now includes both variants. Remote CI completion is not part
of the local validation claim recorded here. Preserve original signed release
snapshots and original payload bytes when making subsequent corrections.
