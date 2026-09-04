# Exact verification bundle for the 12-disk covering theorem

The claimed value is

\[
r_{12}=0.361102963744508644113087702017065180848530565591618515449660087986\ldots.
\]

The exact value is the `r` coordinate of the unique real zero in the rational
box stored in `proof_bundle/cover12_symmetric_kkt_cert.json`.

Run:

```bash
PYTHON_BIN="$(command -v python3)" ./run_verification.sh
```

The mathematical proof is `docs/cover12_proof_zh.md`. The proof is
computer-assisted and has not yet undergone independent peer review.

An independent local replay on September 4, 2026 completed all eight master
verification modules with return code 0. See
`verification/independent_master_replay.log` and `docs/AUDIT_STATUS.md` for the
exact scope and limitations of that check.

Directory layout:

- `proof_bundle/`: exact source, certificates, triangulation data, and module reports;
- `verification/`: master replay logs and the aggregate verification report;
- `docs/`: mathematical proof, candidate image, and audit notes;
- `run_verification.sh`: portable unified entry point.

The final output marker is:

```text
COVER12 EXACT MASTER REPLAY PASSED
```
