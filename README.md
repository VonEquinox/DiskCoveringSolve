# Disk Covering Problem: Computer-Assisted Proofs for n = 11–15

Covering the unit disk with n equal disks, for n = 11, 12, 13, 14, and 15:
global optimality proof claims, exact certificates, and reproducible
verification code. The disk covering problem (also spelled disc covering
problem) asks for the smallest common radius that covers the entire unit disk.

These packages address global optimality, not just numerical local-search
results. The claims remain subject to independent external mathematical and
software review; they do not settle the problem for every n.

Browse the proof packages and verification instructions:
[11 disks](cover11/README.md), [12 disks](cover12/README.md),
[13 disks](cover13/README.md), [14 disks](cover14/README.md),
[15 disks](cover15/README.md).

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

Each exact value is defined by the isolated algebraic root in its rational
Krawczyk certificate, not by the displayed decimal expansion.

## Repository layout

- `cover11/` — eleven-disk proof bundle, documentation, verification records, and runner;
- `cover12/` — twelve-disk proof bundle, documentation, verification records, and runner;
- `cover13/` — thirteen-disk proof bundle, documentation, verification records, and runner;
- `cover14/` — fourteen-disk proof bundle, documentation, verification records, and runner;
- `cover15/` — fifteen-disk primary and alternative bundles, Chinese proofs, and exact replay records;
- `docs/DEVELOPMENT.md` — shared development and verification workflow;
- `requirements.txt` — shared Python dependencies;
- `Makefile` — unified cover11 through cover15 commands.

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
