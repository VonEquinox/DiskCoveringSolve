# Development workflow

This repository contains separate verifier chains for the eleven-disk and
twelve-disk results. Keep their certificates and generated reports isolated:
`cover11/` and `cover12/` are independent proof directories.

## Environment

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Before changing certificates

1. Verify `cover11/proof_bundle/SHA256SUMS` or `cover12/proof_bundle/SHA256SUMS`.
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
- Update the applicable `coverNN/proof_bundle/SHA256SUMS` only after all fixed files are final.

## Adding a new verifier

1. Make it deterministic.
2. Use exact arithmetic for decisive claims.
3. Give every failed assertion a useful context tuple.
4. Add the script to the applicable `coverNN/proof_bundle` verifier chain.
5. Add a machine-readable JSON report when the result links major proof stages.
6. Document the dependency in the applicable `coverNN/docs` directory.
7. Run full replay before merging.
