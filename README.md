# Disk Covering Problem Research

Computer-assisted research and exact verification for covering the unit disk
with congruent disks. The repository currently contains reproducible proof
claims for eleven and twelve disks.

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

Each exact value is defined by the isolated algebraic root in its rational
Krawczyk certificate, not by the displayed decimal expansion.

## Repository layout

- `proof_bundle/` — eleven-disk exact certificates and verifier programs;
- `cover12/` — twelve-disk proof, certificates, verifier programs, and replay logs;
- `docs/METHOD.md` — detailed eleven-disk mathematical method;
- `cover12/cover12_proof_zh.md` — detailed twelve-disk mathematical proof and method;
- `docs/VERIFIER_ARCHITECTURE.md` — eleven-disk verifier dependency graph;
- `docs/DEVELOPMENT.md` — development and verification workflow;
- `verification/` — eleven-disk full replay record.

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
with return code 0. See `cover12/AUDIT_STATUS.md` for the checked scope.

## Integrity

```bash
make integrity
```

The cover12 manifest hashes immutable proof inputs, source files, and fixed
audit documents. Regenerated reports and timing-dependent logs are deliberately
excluded so repeated verification does not invalidate the manifest.

## Status

Both included verifier chains pass locally. These are reproducible
computer-assisted proof claims and should still receive independent external
mathematical and software review.
