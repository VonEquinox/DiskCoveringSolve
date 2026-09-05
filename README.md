# Disk Covering Problem Research

Computer-assisted research and exact verification for covering the unit disk
with congruent disks. The repository currently contains reproducible proof
claims for eleven, twelve, thirteen, and fourteen disks.

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

Each exact value is defined by the isolated algebraic root in its rational
Krawczyk certificate, not by the displayed decimal expansion.

## Repository layout

- `cover11/` — eleven-disk proof bundle, documentation, verification records, and runner;
- `cover12/` — twelve-disk proof bundle, documentation, verification records, and runner;
- `cover13/` — thirteen-disk proof bundle, documentation, verification records, and runner;
- `cover14/` — fourteen-disk proof bundle, documentation, verification records, and runner;
- `docs/DEVELOPMENT.md` — shared development and verification workflow;
- `requirements.txt` — shared Python dependencies;
- `Makefile` — unified cover11 through cover14 commands.

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

## Integrity

```bash
make integrity
```

The cover12, cover13, and cover14 manifests hash immutable proof inputs, source
files, and fixed audit documents. Regenerated reports and timing-dependent logs
are deliberately excluded so repeated verification does not invalidate the
manifest.

## Status

All four included verifier chains pass locally. These are reproducible
computer-assisted proof claims and should still receive independent external
mathematical and software review.
