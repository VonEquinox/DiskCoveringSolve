# Fresh verification record

- Date: 2026-09-04
- Integrity target: 168 fixed files in `proof_bundle/SHA256SUMS`.
- Portable quick audit: passed after the verifier-linkage revisions.
- Full leaf-by-leaf replay: passed with return code 0.
- Terminal marker: `COVER11 FULL LEAF-BY-LEAF REPLAY PASSED`.

The strengthened replay additionally checks Brown's labeled triangulation
counts, the candidate/wheel residual-orbit linkage, exact recomputation of all
six local Kron rows from `weight_nums`, the ninth-gap and 9D-to-8D local-ball
linkage, and the topology of the 29-face upper-bound disk complex.

`full_replay.log` is the complete fresh console log. It was generated from a
temporary copy of `proof_bundle`; the portable wrapper changed only hard-coded
`/mnt/data` path literals and the Python executable path in that copy. The
packaged certificate data and verification algorithms were then executed without
further changes.
