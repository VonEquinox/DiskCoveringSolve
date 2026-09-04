# Exact verification bundle for the 12-disk covering theorem

The claimed value is

\[
r_{12}=0.361102963744508644113087702017065180848530565591618515449660087986\ldots.
\]

The exact value is the `r` coordinate of the unique real zero in the rational box stored in `cover12_symmetric_kkt_cert.json`.

Run:

```bash
PYTHON_BIN="$(command -v python3)" bash verify_all.sh
```

The mathematical proof is `cover12_proof_zh.md`. The proof is computer-assisted and has not yet undergone independent peer review.

An independent local replay on September 4, 2026 completed all eight master
verification modules with return code 0. See `independent_master_replay.log`
and `AUDIT_STATUS.md` for the exact scope and limitations of that check.

The final output marker is:

```text
COVER12 EXACT MASTER REPLAY PASSED
```
