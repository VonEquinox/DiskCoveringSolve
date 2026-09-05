# Development workflow

This repository contains separate verifier chains for eleven through fifteen
disks. Keep their certificates and generated reports isolated: `cover11/`,
`cover12/`, `cover13/`, `cover14/`, and `cover15/` are separate proof directories. Cover15 contains
primary and alternative packages with their own labels, local constants and
enumeration formats.

## Environment

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Before changing certificates

1. Verify the applicable immutable input manifest (`cover11/proof_bundle/SHA256SUMS`, `cover12/proof_bundle/SHA256SUMS`, `cover13/proof_bundle/MANIFEST_SHA256.json`, `cover14/proof_bundle/SHA256SUMS`, or `cover15/INPUT_MANIFEST.json`).
2. Run quick mode and record the baseline output.
3. Identify whether the change is a verifier-only change, a certificate-format
   change, or a mathematical-certificate change.
4. Never replace an exact check with a floating-point tolerance check.

## Required checks for a pull request

Run the cover12 exact master replay whenever any file below `cover12/` changes:

```bash
make cover12-full PYTHON_BIN="$(command -v python)"
```

The required final marker is `COVER12 EXACT MASTER REPLAY PASSED`.

Run the cover13 master replay whenever any file below `cover13/` changes:

```bash
make cover13-full PYTHON_BIN="$(command -v python3)"
```

The required final marker is `COVER13 EXACT MASTER REPLAY PASSED`.  This replay
also requires a C++17 compiler and Boost.Multiprecision headers.

Run the cover14 master replay whenever any file below `cover14/` changes:

```bash
make cover14-full PYTHON_BIN="$(command -v python3)"
```

The required final marker is `COVER14 EXACT MASTER REPLAY PASSED`. This replay
also requires a C++17 compiler and Boost.Multiprecision headers.

When changing the alternative Cover14 bundle or the cross-bundle root linkage,
also run `make cover14-alternative-full`. This uses a temporary working copy so
the alternative archive's 64-file manifest, including submitted reports,
remains valid after replay. The final marker is
`COVER14 ALTERNATIVE MASTER REPLAY AND ROOT LINKAGE PASSED`.

For Cover15 changes, run `make cover15-full`. Both archived packages are copied
to temporary directories, their full censuses are regenerated, and their root
correspondence is checked. To work on one implementation, use
`make cover15-primary-full` or `make cover15-alternative-full`. Check archive
integrity with `make cover15-integrity`; the source-derived
`cover15/INPUT_MANIFEST.json` also covers original submitted reports.

```bash
make integrity
make quick
```

If branch data, stresses, conductances, triangulations, or verifier logic
changes, also run:

```bash
make full
```

## Certificate design rules

- Store integers or exact rational numerator/denominator pairs.
- Recompute derived Kron conductances from original graph weights.
- Check every topology index by an explicit face-set isomorphism witness.
- Separate proposal generation from proof verification.
- State whether a printed decimal is certified, rounded, or diagnostic.
- Keep quick and full success markers distinct.
- Update the applicable immutable manifest only after all fixed files are final.

## Adding a new verifier

1. Make it deterministic.
2. Use exact arithmetic for decisive claims.
3. Give every failed assertion a useful context tuple.
4. Add the script to the applicable `coverNN/proof_bundle` verifier chain.
5. Add a machine-readable JSON report when the result links major proof stages.
6. Document the dependency in the applicable `coverNN/docs` directory.
7. Run full replay before merging.
