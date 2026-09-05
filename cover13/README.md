# Exact verification bundle for the 13-disk covering theorem

The claimed optimum is the algebraic number

\[
r_{13}=0.346645456927389643467869738193878789537187452966783984053969\ldots.
\]

More precisely, if \(t_*\) is the radius-squared coordinate of the unique real
zero in the rational box stored in
`proof_bundle/core/cover13_krawczyk_cert.json`, then the theorem claims
\(r_{13}=\sqrt{t_*}\).  The verified strict interval is

\[
0.346645456927389643467869738193878789537187452966783984053969
<r_{13}<
0.346645456927389643467869738193878789537187452966783984053970.
\]

## Replay

Requirements:

- Python 3.10 or later;
- a C++17 compiler (`g++` or `clang++`);
- Boost.Multiprecision headers.

On macOS with Homebrew:

```bash
brew install boost
PYTHON_BIN="$(command -v python3)" ./run_verification.sh
```

On Debian/Ubuntu:

```bash
sudo apt-get install g++ libboost-dev
PYTHON_BIN="$(command -v python3)" ./run_verification.sh
```

The final marker is:

```text
COVER13 EXACT MASTER REPLAY PASSED
```

The runner first checks all 28 immutable input hashes, then reruns every Python
checker, recompiles the C++ orbit auditor, and independently audits all four
triangulation families.  Generated reports are written to
`proof_bundle/reports/` and are ignored by Git.

## Contents

- `docs/cover13_proof_zh.md` - detailed Chinese mathematical proof;
- `proof_bundle/` - exact certificates, verifier source, enumeration data, and
  the fail-closed master runner;
- `verification/` - submitted reports and the independent replay performed on
  September 5, 2026;
- `run_verification.sh` - portable repository entry point;
- `verify_manifest.py` - independent SHA-256 manifest check.

## Verified scope

The replay checks the 119-variable rational Krawczyk root certificate, the
strict upper construction, cross-file interfaces, local anchor isolation, the
candidate global branch tree, the exact metric partition, 1,760 noncandidate
residual certificates, negative fail-closed tests, and orbit completeness for
all 308,198 symmetry classes.

This is a replayable computer-assisted proof claim.  It is not a Lean/Coq
formalization and has not yet undergone independent peer review.
