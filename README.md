# DiskCoveringSolve

Covering the unit disk with n equal disks, for n = 11 through 20:
global optimality proof claims, exact certificates, and reproducible
verification code. The disk covering problem (also spelled disc covering
problem) asks for the smallest common radius that covers the entire unit disk.

These packages address global optimality, not just numerical local-search
results. The claims remain subject to independent external mathematical and
software review; they do not settle the problem for every n.

Browse the proof packages and verification instructions:
[11 disks](cover11/README.md), [12 disks](cover12/README.md),
[13 disks](cover13/README.md), [14 disks](cover14/README.md),
[15 disks](cover15/README.md), [16 disks](cover16/README.md),
[17 disks](cover17/README.md), [18 disks](cover18/README.md),
[19 disks](cover19/README.md), [20 disks](cover20/README.md).

The repository also contains a [Chinese proof of a single-exponential exact
algorithm for arbitrary disk counts](theory/single_exponential_exact_algorithm_zh.md).
This is a theoretical algorithm and bit-complexity result, not an implemented
end-to-end solver or a table of computed optimal radii for every n. It is
independent of the fixed-n certificate packages and remains open to external
mathematical review.

## Claimed results

For

\[
r_n=\inf_{c_1,\ldots,c_n\in\mathbb R^2}
\max_{\lVert x\rVert\le1}\min_i\lVert x-c_i\rVert,
\]

the included certificate packages claim:

| Disks | Certified value |
|---|---|
| 11 | \(r_{11}=0.37998385311983868972226160613609321430871895\ldots\) |
| 12 | \(r_{12}=0.36110296374450864411308770201706518084853056\ldots\) |
| 13 | \(r_{13}=0.34664545692738964346786973819387878953718745\ldots\) |
| 14 | \(r_{14}=0.33173203427623412255781554880240953827094275\ldots\) |
| 15 | \(r_{15}=0.31814293085926283635949118164993500518666386\ldots\) |
| 16 | \(r_{16}=0.30821980189861784887010326922181483200610642\ldots\) |
| 17 | \(r_{17}=0.29853158925642805986984377442398070373438760\ldots\) |
| 18 | \(r_{18}=0.29016771764056672911315452560106171830628088\ldots\) |
| 19 | \(r_{19}=1/\sqrt{13}=0.27735009811261456100917086672849968817317665\ldots\) |
| 20 | \(r_{20}=0.27084817843997736519402387625893583587132712\ldots\) |

For n = 11 through 18, the exact value is defined by the isolated algebraic
root in its rational Krawczyk certificate, not by the displayed decimal
expansion. Cover19 instead uses the explicit radical \(1/\sqrt{13}\) and
exact geometry in \(\mathbb Q(\sqrt{3})\); it needs no root-isolation certificate.
Cover20 defines its exact value by a 198-variable integer quadratic system
and two nested rational contraction boxes.

## Repository layout

- `cover11/` — eleven-disk proof bundle, documentation, verification records, and runner;
- `cover12/` — twelve-disk proof bundle, documentation, verification records, and runner;
- `cover13/` — thirteen-disk proof bundle, documentation, verification records, and runner;
- `cover14/` — fourteen-disk proof bundle, documentation, verification records, and runner;
- `cover15/` — fifteen-disk primary and alternative bundles, Chinese proofs, and exact replay records;
- `cover16/` — sixteen-disk proof, exact certificates, census source, and full replay records;
- `cover17/` — seventeen-disk proof, exhaustive peeling sources, exact certificates, and full replay records;
- `cover18/` - eighteen-disk primary and alternative proofs, exact root linkage, discovery sources, and full replay records;
- `cover19/` - nineteen-disk primary and alternative proofs, exact geometry and stress linkage, discovery sources, and full replay records;
- `cover20/` - twenty-disk ZIP proof, exact certificates with lossless large-file parts, discovery sources, and fresh replay records;
- `theory/` - Markdown proofs of general algorithms, separate from fixed-n certificate replays;
- `docs/DEVELOPMENT.md` — shared development and verification workflow;
- `requirements.txt` — shared Python dependencies;
- `Makefile` — unified cover11 through cover20 commands.

Manuscript PDF/LaTeX files are intentionally excluded. This repository is for
proof development, exact certificates, and reproducible verification.

## Environment

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

## Eleven-disk verification

Quick structural audit:

```bash
make cover11-quick PYTHON_BIN="$(command -v python)"
```

Full leaf-by-leaf replay:

```bash
make cover11-full PYTHON_BIN="$(command -v python)"
```

Expected full marker:

```text
COVER11 FULL LEAF-BY-LEAF REPLAY PASSED
```

## Twelve-disk verification

The twelve-disk verifier checks the rational Krawczyk box, symbolic KKT
identities, upper cover, six exhaustive triangulation families, exact Farkas
certificates, the complete residual partition, 131 graph-energy certificates,
and the candidate global-convexity certificate.

```bash
make cover12-full PYTHON_BIN="$(command -v python)"
```

Expected marker:

```text
COVER12 EXACT MASTER REPLAY PASSED
```

A fresh local replay on September 4, 2026 completed all eight cover12 modules
with return code 0. See `cover12/docs/AUDIT_STATUS.md` for the checked scope.

## Thirteen-disk verification

The thirteen-disk master verifier checks the 119-variable root certificate,
the exact upper construction, the candidate and noncandidate branch trees,
the metric partition, cross-stage interfaces, fail-closed negative tests, and
all four orbit-enumeration families.  It requires a C++17 compiler and
Boost.Multiprecision headers in addition to Python.

```bash
make cover13-full PYTHON_BIN="$(command -v python3)"
```

Expected marker:

```text
COVER13 EXACT MASTER REPLAY PASSED
```

A fresh local replay on September 5, 2026 passed all stages and audited all
308,198 topology orbits. See `cover13/docs/AUDIT_STATUS.md`.

## Fourteen-disk verification

The fourteen-disk master verifier checks the 116-variable root certificate,
exact symmetry and upper-cover identities, anchor isolation, the candidate and
noncandidate branch trees, the exact metric partition, fail-closed negative
tests, and all four orbit-enumeration families. It requires a C++17 compiler
and Boost.Multiprecision headers in addition to Python.

```bash
make cover14-full PYTHON_BIN="$(command -v python3)"
```

Expected marker:

```text
COVER14 EXACT MASTER REPLAY PASSED
```

A fresh local replay on September 5, 2026 passed all stages and audited all
1,313,024 topology orbits. See `cover14/docs/AUDIT_STATUS.md`.

An alternative Cover14 certificate package also passes full replay, with smaller
branch trees and a larger anchor-isolation domain. An exact root-box transport
identifies its candidate with the primary package. Run
`make cover14-alternative-full`; see `cover14/docs/ALTERNATIVE_AUDIT.md` for the
comparison and shared dependencies.

## Fifteen-disk verification

Cover15 includes two related certificate chains for the same algebraic
candidate. Each regenerates and independently audits all 11,950,884 topology
orbits, verifies the 148-variable root and upper cover, and checks all 43,013
noncandidate cases and the entire candidate domain. A separate exact root-box
transport binds their different coordinate labels.

```bash
make cover15-full PYTHON_BIN="$(command -v python3)"
```

The final marker is `COVER15 REQUESTED REPLAYS AND ROOT LINKAGE PASSED`.
Individual chains are available through `cover15-primary-full` and
`cover15-alternative-full`. Allow several GB of disk and memory and several
minutes for regeneration. See `cover15/docs/AUDIT_STATUS.md` for provenance,
the comparison, replay results and review limits.

## Sixteen-disk verification

Cover16 regenerates all 53,059,205 topology orbits and audits them with a
separate breadth-first normalization and exact labelled counting recurrence.
The 140-variable root, 42-triangle upper cover, local matrix, four- and six-rod
angle constraints, 175 noncandidate types and full candidate angle domain
are checked using exact arithmetic.

```bash
make cover16-full PYTHON_BIN="$(command -v python3)"
```

Use Python 3.10+ and a C++17 compiler with Boost headers. The runner uses one
census worker by default and discards temporary enumeration tables afterward.
The final marker is `COVER16 EXACT MASTER REPLAY PASSED`.
See [Cover16](cover16/README.md) and its [audit record](cover16/docs/AUDIT_STATUS.md)
for the proof, resource requirements and fresh replay evidence.

## Seventeen-disk verification

Cover17 uses exhaustive root-face peeling with necessary-condition pruning
on partial triangulations. Two implementations regenerate and compare the
complete surviving sets for all six families. There are 1,344 surviving
orbits, not 1,344 total unpruned triangulations. All 1,343 noncandidate types
and the entire candidate angle domain are then checked exactly.

```bash
make cover17-full PYTHON_BIN="$(command -v python3)"
```

The final marker is `COVER17 EXACT MASTER REPLAY PASSED`. This chain needs
Python 3.10+ and a C++17 compiler, but no Boost or third-party Python packages.
See [Cover17](cover17/README.md) and the [audit record](cover17/docs/AUDIT_STATUS.md)
for the 157-variable root, 44-triangle upper cover, pruning-completeness
argument and fresh replay results.

## Eighteen-disk verification

Cover18 includes two related certificate chains. Each checks the 174-variable
root, the 46-triangle upper cover, exhaustive peeling over seven families,
all 16,219 noncandidate types and the full candidate angle domain. The 16,220
surviving orbits are not the total unpruned topology count. An additional
exact root transport and full survivor-set comparison link the two archives.

```bash
make cover18-full PYTHON_BIN="$(command -v python3)"
```

Python 3.10+ and a C++17 compiler suffice; no third-party Python dependency or
Boost is needed. The final marker is
`COVER18 REQUESTED REPLAYS AND EXACT LINKAGE PASSED`.
The alternative package also includes optional numerical discovery sources,
which are not used for acceptance. See [Cover18](cover18/README.md) and its
[audit record](cover18/docs/AUDIT_STATUS.md) for both proofs and actual results.

## Nineteen-disk verification

Cover19 includes two related exact certificate chains for
\(r_{19}=1/\sqrt{13}\). Both check the upper construction, exhaustive peeling
over seven families, and all 24,126 noncandidate types. The primary handles
the candidate graph by global energy convexity; the alternative uses exact
anchor isolation and a complete candidate branch tree. The 24,127 surviving
orbits are survivors after necessary-condition pruning, not an unpruned census.

```bash
make cover19-full PYTHON_BIN="$(command -v python3)"
```

Python 3.10+ and a C++17 compiler suffice, without third-party Python packages
or Boost. The final marker is
`COVER19 REQUESTED REPLAYS AND EXACT LINKAGE PASSED`.
An additional checker compares both exact constructions, their 42-rod stresses,
and every surviving orbit. See [Cover19](cover19/README.md) and its
[audit record](cover19/docs/AUDIT_STATUS.md) for both Chinese proofs, the
degenerate upper-cover faces, the omitted central stress vertex, and fresh
full replay evidence.

## Twenty-disk verification

Cover20 checks the 198-variable isolated root, 51-triangle upper cover,
102-coordinate anchor-isolation matrix, two exhaustive peeling searches over
eight families, all 562,168 noncandidate cases and the entire candidate tree.
The 562,169 orbit classes are survivors after safe partial-state pruning,
not the unpruned topology total.

```bash
make cover20-full PYTHON_BIN="$(command -v python3)"
```

The runner reconstructs one large certificate from three lossless parts,
checks its original ZIP hash, and replays the unchanged original master in a
temporary directory. No external downloads or Git LFS are required.
Python 3.10+, a C++17 compiler, and several GB of available disk/memory suffice.
The final marker is `COVER20 EXACT MASTER REPLAY AND INPUT RECONSTRUCTION PASSED`.
See [Cover20](cover20/README.md) and its [audit record](cover20/docs/AUDIT_STATUS.md).
Only the ZIP's matching Chinese proof is included, not a separately supplied
Markdown version with different certificate definitions.

## Integrity

```bash
make integrity
```

The cover12, cover13, and cover14 manifests hash immutable proof inputs, source
files, and fixed audit documents. Regenerated reports and timing-dependent logs
are deliberately excluded so repeated verification does not invalidate the
manifest.

Cover15 preserves both original archives' reports and checks all archived
files. Its replay runs in disposable copies so regenerated reports and census
files do not change those archived hashes.

Cover16 likewise preserves all 90 archived files, including submitted reports.
Its complete census replay runs in a disposable copy, and fresh audit evidence
is stored separately from the original archive payload.

Cover17 preserves all 46 archived files, including supplied stage reports.
Its runner removes stale reports only in a disposable copy before executing
every accepting stage. Fresh evidence is stored separately in `verification/`.

Cover18 preserves all 107 original files across its two archives. Its runner
removes old reports only in disposable copies, checks every requested stage,
and stores new evidence separately from the submitted payloads.

Cover19 follows the same policy for all 89 original files. Its two standalone
Chinese proofs are also checked against the originals. Supplied reports are
not counted as freshly executed checks; each full run regenerates its reports
in disposable copies and then checks exact cross-bundle linkage.

Cover20 preserves all 51 original logical files. Fifty are stored directly;
the largest compressed forest is stored as three byte ranges in `proof_parts/`.
Integrity verification checks each part's size and hash and the concatenated
file's original hash. The full runner restores the exact original file only
in a disposable copy before executing any proof checks.

## Status

The included results are reproducible
computer-assisted proof claims and should still receive independent external
mathematical and software review.

## Model provenance

As reported by the repository maintainer, the proof packages were produced
through the following model sessions (model names are recorded as supplied):

| Problems | Model / interface | Reported thinking time |
|---|---|---|
| Cover11, Cover12 | GPT 5.6 Pro, web interface | Extended thinking sessions; exact durations not recorded here |
| Cover13, Cover14, Cover15 | GPT 6 astra Pro | No more than 50 minutes per problem |

These are maintainer-reported model thinking times, not verifier runtimes or
total project durations. Model provenance is not evidence of correctness;
the mathematical arguments and reproducible certificate checks are the basis
for assessing the proof claims.
