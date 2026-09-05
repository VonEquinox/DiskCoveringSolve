# Cover17 audit record

## Input provenance

Reviewed on September 5, 2026:

| Input | SHA-256 |
|---|---|
| `r17_verified_proof.zip` | `58336cd10c181d93d67ebdd9e0d20de3d7522b8402a6e82cde6840533778a357` |
| `r17_proof_zh.md` | `fb429e7b489fdda94be5c5a099a0c950a23b02f78a09b5bbad7fa24dfd367b92` |

Archive paths, entry uniqueness, CRCs and all 30 accepting-input hashes in
the supplied master report passed before execution. The standalone Chinese
proof is byte-identical to `PROOF_zh.md`. This archive has no separate
SHA256SUMS file; the repository adds `INPUT_MANIFEST.json` covering all
46 original files, including supplied stage reports. None of the original
proof text, certificate data or accepting source files was edited.

The root certificate SHA-256 is
`f520c34816aab807fd115933a1159151ec7e71e4b8fb3b9b9d85d7c554050ac2`.

## Fresh full replays

Both runs below recompiled the two C++ enumeration implementations, regenerated
all six families and replayed all accepting stages. Neither used the supplied
master's success flag as acceptance evidence.

| Run | Command | Result | Master elapsed seconds |
|---|---|---|---:|
| Fresh archive extraction | `python3 -S -B verify_all.py` | PASS | 13.985297918319702 |
| New repository entry point | `python3 -S -B cover17/run_verification.py --report-dir cover17/verification/repository_replay` | PASS | 13.674449920654297 |

The local interpreter was Homebrew Python 3.14.6, with a C++17 compiler.
These are observed local timings, not controlled benchmarks or model thinking
times. No Boost or external Python packages were required. Both runs returned
code 0, verified 1,344 surviving orbits and reported zero unresolved leaves.

Fresh reports and logs are in `../verification/independent_MASTER_VERIFIED.json`,
`independent_full_replay.log`, `independent_reports/`, `repository_full_replay.log`
and `repository_replay/`. The exported repository replay includes all freshly
generated core stage summaries. Original submitted reports remain separate
inside `proof_bundle/`.

## What the enumeration count means

This is not the Cover16-style unpruned orbit census. The new exhaustive
algorithm prunes partial triangulations using necessary angle inequalities.
Its written completeness proof must justify both root-face choices and the
safety of every pruning condition.

| Family | Primary partial states | Independent partial states | Rooted surviving records | Surviving orbits |
|---|---:|---:|---:|---:|
| B11 I6 | 161,886 | 102,179 | 10,439 | 489 |
| B12 I5 | 177,422 | 242,611 | 14,583 | 631 |
| B13 I4 | 90,653 | 360,159 | 4,992 | 200 |
| B14 I3 | 44,600 | 477,303 | 0 | 0 |
| B15 I2 | 26,820 | 429,475 | 0 | 0 |
| B16 I1 | 15,293 | 204,649 | 512 | 24 |
| Total | 516,674 | 1,816,376 | 30,526 | 1,344 |

The exact surviving-set partition is 1,343 noncandidate types plus one
candidate. The candidate tree contains 2,787 nodes: 1,393 splits, 1,357 force
leaves and 37 local leaves. The noncandidate forest has 4,075 nodes:
1,366 splits and 2,709 force leaves. No branch is left unresolved.

## Claim and evidence alignment

| Claim | Evidence reviewed and replayed | Status |
|---|---|---|
| Precisely defined radius | Rebuilt 157-variable integer KKT system, two nested contraction boxes and exact decimal comparisons | PASS |
| Feasible full unit-disk cover | 44 oriented triangles, 12 caps, 13 exact active faces, 7 strictly slack faces and boundary winding argument | Finite checks passed; written argument reviewed |
| Arbitrary smaller covers enter six families | Strict-buffer ordinary Voronoi reduction, irreducibility/private points, no-return, simple dual and stellar padding | Written argument reviewed; not proof-assistant formalized |
| Partial-state pruning cannot lose a counterexample | Root-face state invariant/induction and necessary short-path difference inequalities | Written argument and both implementations reviewed |
| Complete surviving sets agree | Fresh stack/DFS/incremental-closure search versus FIFO/BFS/full-Floyd search, with exact full-code set comparison | PASS for all six families |
| All residual continuous domains are excluded | Full graph/isomorphism linkage, integer forces, complete tree reachability and anchor-only isolation | PASS |
| Rejection behavior and unpruned enumeration are tested | 30 supplied negative tests, 22 supplied unpruned counts, 6 additional sanitized unpruned cases | PASS; finite testing is not a universal software guarantee |

## Mathematical and source review

The full Chinese proof and all acceptance-source modules were read. Its
sequence is coherent: algebraic definition and upper construction, arbitrary
cover reduction, safe partial-state inequalities, exhaustive peeling,
survivor linkage, global force leaves and certified local leaves. It does
not assume candidate symmetry, twelve boundary disks for an arbitrary cover,
or optimality at any smaller n. It does not claim to classify all optimizers.

The geometric reduction uses ordinary Voronoi cells after a strict-radius
buffer and an irreducible subcover. The no-return proof uses a private point
and a convex circular segment. The dual and stellar-padding arguments lead
to B=11 through 16 with I=17-B; inserted coincident auxiliary centers are
allowed by the later necessary conditions.

For peeling, installing a root face either splits a region along a boundary
third vertex or exposes one fresh interior vertex. Internal-vertex budgets
are exhausted over all possible splits. Already installed nonboundary edges
are rejected when their reuse would create parallel edges. The remaining
face count decreases by one at each step. A genuine completion therefore
follows one branch, up to deterministic relabeling of interior vertices.
Only edges already forced in every completion are used for pruning.

The primary implementation maintains a closed integer difference-bound
matrix incrementally. The second reconstructs adjacency from installed faces
and runs full Floyd-Warshall with wider integer storage. Both disambiguate
short/long angular branches before imposing directed bounds. The audit's
final comparison is an exact set comparison after independent BFS coding,
not just an equality of totals. The written induction is still essential;
two related implementations agreeing is not by itself a completeness proof.

The upper complex has (V,E,F)=(29,72,44). Its positive triangle orientations,
opposite internal-edge incidences and prescribed boundary cycle give the
winding-number coverage argument. The root-linked angle intervals explicitly
bound each consecutive gap, including the wraparound gap, below a quarter
turn before the circular-cap argument is applied.

There are 82 free geometry coordinates and 74 geometric constraint rows.
The positive-weight graph has 63 rods and embeds explicitly into the full
49-node, 84-rod candidate graph. The local isolation radius is 27/250, with
exact margin 499/78125000. Its premise controls anchors only; it does not
assume the other free nodes are close to the candidate. Candidate orbit
B12 I5 index 630, combined residual index 1119, is linked by the checked
center permutation `[11,0,1,2,3,4,5,6,7,8,9,10,16,15,13,12,14]`.

The exact trigonometric source is identical to the Cover16 implementation.
Other interval/root/force components also follow related earlier patterns.
These are related proof packages, not independent teams or third-party peer
reviews. No specific mathematical or certificate error blocking publication
was identified during this review.

## Additional checks

The supplied no-pruning regressions only exercise at most three interior
vertices. `test_peeling_counts.py` supplements them with the following six
cases, using both C++ implementations, no angle pruning and the exact labelled
count recurrence. Both executables were compiled with AddressSanitizer and
UndefinedBehaviorSanitizer; all tests passed without sanitizer diagnostics.

| B | I | Rooted count | Labelled count |
|---|---|---:|---:|
| 3 | 4 | 68 | 1,632 |
| 3 | 5 | 399 | 47,880 |
| 3 | 6 | 2,530 | 1,821,600 |
| 4 | 4 | 570 | 13,680 |
| 4 | 5 | 3,542 | 425,040 |
| 4 | 6 | 23,400 | 16,848,000 |

Evidence: `../verification/extra_unpruned_counts_verified.json` and
`extra_unpruned_counts.log`. These supplement the proof; they do not replace
the large-n completeness induction. Eight synthetic runner tests additionally
check input hashes, missing files, proof mismatch, manifest paths, temporary
copy/export cleanup, skipped accepting stages, stale reports and protected
output locations. These are not counted as further mathematical full replays.

After integration, Cover11 through Cover17 passed `make integrity`, all eight
runner tests passed, and a byte-for-byte comparison confirmed that all 46
archived files match the original ZIP. The checks are recorded in
`../verification/repository_integrity_and_runner_tests.log`.

The trust boundary remains the written geometry, topology and peeling
induction, Python exact arithmetic, and the C++ toolchain. External peer
review, an independently developed external verifier and Lean/Coq
formalization have not been completed. Local checks do not assert a remote
GitHub Actions run has already passed.

The proof cites [Erich Friedman's covering table](https://erich-friedman.github.io/packing/circovcir/)
for historical construction attribution to Jeremy Tan (2018), not as a
premise for optimality. No first-discovery claim for that numerical candidate
is made here. No arXiv manuscript files are included.
