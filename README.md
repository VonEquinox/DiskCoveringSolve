# Disk Covering Problem Research

Computer-assisted research and exact verification for covering the unit disk
with eleven congruent disks.

The repository contains the complete verification code, rational
certificates, triangulation data, branch trees, replay logs, and detailed
method documentation. Manuscript PDF/LaTeX files are intentionally excluded;
this is the development and verification repository.

## Claimed result

Let

\[
r_{11}=\inf_{c_1,\ldots,c_{11}\in\mathbb R^2}
\max_{\|x\|\le1}\min_i\|x-c_i\|.
\]

The certificate package claims

\[
r_{11}=0.37998385311983868972226160613609321430871895\ldots.
\]

The exact value is not defined by the decimal. It is
`r_* = sqrt(t_*)`, where `t_*` is the `t`-coordinate of the unique KKT
root inside the rational box certified by
`proof_bundle/cover11_full_kkt_certificate.json`.

## Proof strategy

The complete method is described in [`docs/METHOD.md`](docs/METHOD.md). In
outline:

1. isolate an exact algebraic KKT candidate using a rational Krawczyk test;
2. certify that the candidate covers the unit disk;
3. reduce any hypothetical smaller cover to an ordinary restricted Voronoi
   diagram in general position;
4. prove that its dual is one of 3843 simple disk-triangulation orbits;
5. eliminate 3788 orbits by exact Farkas certificates;
6. eliminate the wheel and 53 ordinary residuals by rational graph energies;
7. handle the candidate topology by global branch-and-bound plus an
   eight-dimensional strict local-minimum certificate.

## Repository layout

- `proof_bundle/` - exact proof data and verification programs;
- `verification/` - full replay log and machine-readable result summary;
- `docs/METHOD.md` - detailed mathematical and computational method;
- `docs/VERIFIER_ARCHITECTURE.md` - certificate dependency graph and trusted
  computational base;
- `docs/DEVELOPMENT.md` - development workflow and extension guide;
- `run_verification_portable.sh` - unified portable replay entry point;
- `requirements.txt` - Python dependencies.

## Integrity

```bash
cd proof_bundle
sha256sum -c SHA256SUMS
```

The proof-bundle manifest contains 168 fixed files.

## Quick audit

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
./run_verification_portable.sh quick
```

The final marker must be:

```text
COVER11 QUICK AUDIT PASSED (STORED LEAF SUMMARIES ONLY)
```

Quick mode checks the central exact certificates, enumeration, Brown counts,
Farkas structure, topology linkage, local Kron recomputation, the 9-gap/8D
link, the upper complex, and stored leaf summaries. It does not recompute every
branch leaf.

## Full leaf-by-leaf replay

```bash
./run_verification_portable.sh full
```

The final marker must be:

```text
COVER11 FULL LEAF-BY-LEAF REPLAY PASSED
```

See `verification/full_replay.log` for the latest complete replay transcript.

## Current verification totals

| Stage | Result |
|---|---:|
| Triangulation orbits | 3843 |
| Farkas-eliminated orbits | 3788 |
| Candidate direct leaves | 100834 |
| Candidate refined external leaves | 153052 |
| Candidate local leaves | 1953 |
| Ordinary residual topologies | 53 |
| Ordinary residual strict leaves | 13789 |

## Status

The included strengthened verifier has completed a fresh full replay with
return code 0. The result is computer-assisted and should receive independent
external mathematical and software audit.
