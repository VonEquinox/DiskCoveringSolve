# Cover16 audit record

## Inputs and scope

Reviewed on September 5, 2026. The accepted archive and separate manuscript:

| Input | SHA-256 |
|---|---|
| `r16_verified_proof.zip` | `fc6f8ddd61ad14cba29c6fa732c0eb3727dbf644fbe32e9c23bfecdc1321003e` |
| `r16_proof_zh.md` | `532c3d5513967a6a0f682e596212f60a6231aa0e564031a5912bf9e0fb70f3d9` |

ZIP path safety, entry uniqueness, CRC checks and the supplied checksum list
passed before execution. The separate Chinese manuscript is byte-identical
to `PROOF_zh.md`. All 90 original archive files, including discovery programs
and submitted reports, are preserved unchanged in `../proof_bundle/`.
The repository's `INPUT_MANIFEST.json` hashes every archived file.

A second file named `r16_verified_proof(1).zip` was also submitted during this
review. At inspection it was 13,764,650 bytes, with SHA-256
`13868548b43b6ad512e2f90346473e53a57446f98a0bb1a3d0c9b0c7b86e42c6`.
It had a ZIP local header but lacked the end-of-central-directory record and
could not be extracted. It is not included as a verified alternative. This
is a file-integrity problem, not a finding about an alternative mathematical
proof; no comparison of its certificates is claimed.

## Fresh full replay

The accepted archive was extracted into an empty working directory, without
any census tables. The following command was executed locally using Homebrew
Python 3.14.6, the C++17 toolchain and Homebrew Boost headers:

```sh
CPLUS_INCLUDE_PATH="$(brew --prefix boost)/include" \
  /opt/homebrew/bin/python3 -S -B verify_all.py --jobs 2
```

It returned code 0 and printed `R16 COMPLETE CERTIFICATE CHAIN VERIFIED`.
The fresh master reports `regenerated_and_independently_audited`, elapsed
time **530.712651014328 seconds**, and **zero unresolved leaves**. This is
a local observed runtime, not a benchmark or model thinking-time measurement.

The supplied report used separately generated tables and a reuse-mode audit;
that report was not substituted for the new from-scratch replay. The fresh
master and top-level log are `../verification/independent_MASTER_VERIFIED.json`
and `../verification/independent_full_replay.log`. Stage reports, generation
summaries, census logs and audited residual records are retained under
`../verification/independent_reports/`.

| Family | Regenerated and audited orbits | Exact labelled mass |
|---|---:|---:|
| B11 I5 | 29,372,227 | 77,519,922,480 |
| B12 I4 | 14,581,757 | 8,394,640,800 |
| B13 I3 | 6,294,651 | 981,608,628 |
| B14 I2 | 2,230,735 | 124,807,200 |
| B15 I1 | 579,835 | 17,383,860 |
| Total | 53,059,205 | |

The partition is exactly 51,464,643 direct metric exclusions plus 1,594,386
Farkas exclusions plus 175 noncandidate residual types plus one candidate.
The candidate tree has 10,347 nodes: 5,173 splits, 5,052 force leaves and
122 local leaves. The noncandidate forest has 713 nodes: 269 splits and
444 force leaves. Every strict acceptance margin is checked rationally.

## Claim and evidence alignment

| Claim | Evidence checked | Status |
|---|---|---|
| A well-defined algebraic radius | Rebuilt 140-variable KKT system; contraction in nested rational boxes; exact decimal interval comparisons | Full replay passed |
| Candidate covers the unit disk | 42 oriented triangles, 12 circular caps, 10 exact active faces and 8 strict slack faces | Full replay passed; winding argument reviewed |
| Every smaller cover enters the finite families | Strict buffered perturbation, irreducible ordinary Voronoi cells, no-return, dual disk and stellar padding in Section 4 | Written argument reviewed, not proof-assistant formalized |
| Finite census has no omitted or duplicate type | Validity, BFS canonicalization, stabilizers, distinct orbit codes and exact rooted-counting mass | Full replay passed |
| Every residual continuous domain is covered | Exact metric/graph linkage, all reachable tree nodes, complete residual set, force inequalities and local isolation | Full replay passed |
| Invalid certificates are rejected in tested cases | 20 supplied adversarial tests | Full replay passed; not exhaustive software testing |

## Mathematical and implementation review

The complete Chinese manuscript and acceptance-source files were read. Its
argument proceeds from an explicitly defined root to a feasible upper cover,
then from arbitrary smaller covers to a complete finite census, exact metric
cuts and exhaustive continuous-domain certificates. No particular symmetry
or twelve-boundary/four-interior structure is assumed for arbitrary covers.
The proof does not depend on optimality results for smaller values of n.

The no-return argument uses a private point of an irreducible cover and a
convex circular segment. It does not assert that the entire connecting arc
is disjoint from the original cell. Genericization uses ordinary Voronoi
cells after a strict-radius buffer, so nearest-center containment follows
directly. Stellar padding keeps the boundary count fixed and permits
coincident auxiliary centers; the later necessary conditions allow this.
The dual-triangulation and rooted peeling/counting arguments remain written
mathematics in the trust boundary.

The upper verifier proves positive orientation, opposite directed internal
edge incidences, one prescribed boundary cycle, Euler characteristic and
vertex-link conditions. The manuscript uses cancellation of oriented triangle
boundaries and winding numbers to conclude coverage; it does not rely only
on total area. Exact common-point and circular-cap inequalities supply the
disk coverage. The upper complex has (V,E,F) = (28,69,42).

The local problem has 74 geometry coordinates, 65 geometric constraints,
54 positive-weight rods, and eleven independent anchor angles. The local
matrix is positive definite after the stated rational perturbation budget.
Its anchor-only radius is 1/10, with isolation margin 43/500000. The proof
does not assume the other coordinates are close. The full candidate graph
has 46 nodes and 78 rods; the verifier checks the exact embedding of the
54-rod activity subgraph. Deleting the other constraints weakens the
minimax problem, so the subgraph lower bound suffices.

The new six-rod constraints distinguish short and long angular branches
before accepting directed arc bounds. The six-rod tag is stripped from the
gap mask and assigned its own right-hand-side constant in both the C++ and
Python Farkas checks. Force divergence is checked as integer identities;
arc support and norm bounds are outward rounded. Floyd-Warshall contraction,
strict interior splits and reachability checks cover the entire initial
angle domain. Residual indices and complete face sets are linked explicitly,
including candidate orbit B12 I4 index 13,870,879, residual index 80.

The census generator uses DFS rooted traversal. The auditor uses BFS rooted
traversal and a different boundary invariant; shared storage/hash semantics
compare full codes on collisions. Completeness comes from validity,
distinctness and exact orbit mass, not from trusting flip-search termination.
This is algorithmic cross-checking inside the supplied proof package, not
an independently developed external verifier.

## Repository checks and limits

The new runner checks immutable inputs, runs the unmodified master on a
temporary copy with regeneration enabled, rejects stale or missing success
reports, optionally exports fresh reports, and discards generated census
tables. Eight lightweight tests cover hash tampering, missing inputs,
manuscript mismatch, unsafe manifest paths, temporary-copy/export cleanup,
reused-census rejection, stale-report rejection and protected output paths.
They use synthetic fixtures and do not count as another mathematical replay.
The full mathematical run above used the original master directly.

After integration, `make integrity` passed for Cover11 through Cover16 and
all eight runner tests passed. The log is
`../verification/repository_integrity_and_runner_tests.log`. A separate
byte-for-byte comparison confirmed that all 90 archived files still match
the original ZIP. No original certificate or acceptance source was edited.

No specific certificate inconsistency or mathematical error blocking
publication was identified in this review. Successful program execution
alone is not a formal proof of the written geometric reduction or a guarantee
against all software defects. External peer review and Lean/Coq formalization
have not been completed. CI is configured separately; local success does
not assert that a remote CI run has already passed.

The manuscript cites [Erich Friedman's Circles Covering Circles table](https://erich-friedman.github.io/packing/circovcir/)
for construction background. The table credits the displayed n=16 construction
to Jeremy Tan (2018). The submission's claimed contribution is global
optimality with exact certificates, not discovery of a previously unknown
numerical candidate. This background table is not an input to acceptance.
