# Exact verification bundle for the 14-disk covering theorem

The claimed optimum is the algebraic number

\[
r_{14}=0.331732034276234122557815548802409538270942754829\ldots.
\]

More precisely, if \(t_*\) is coordinate 62 (zero-based) of the unique real
zero in the rational box stored in `proof_bundle/core/root14.json`, then the
theorem claims \(r_{14}=\sqrt{t_*}\). The verified strict interval is

\[
0.331732034276234122557815548802409538270942754829017478062828
<r_{14}<
0.331732034276234122557815548802409538270942754829017478062829.
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
COVER14 EXACT MASTER REPLAY PASSED
```

The runner first checks all 37 fixed source, certificate, enumeration, and
proof-document hashes. It then reruns every Python checker, recompiles the C++
orbit auditor, and independently audits all four triangulation families.
Generated reports are written to `proof_bundle/reports/` and are ignored by
Git.

## Contents

- `docs/cover14_proof_zh.md` - detailed Chinese mathematical proof;
- `proof_bundle/` - exact certificates, verifier and discovery source,
  enumeration data, and the fail-closed master runner;
- `verification/` - submitted reports and the independent replay performed on
  September 5, 2026;
- `run_verification.sh` - portable repository entry point;
- `verify_manifest.py` - independent SHA-256 manifest check.

## Verified scope

The replay checks the 116-variable rational root certificate, exact symmetry
identities, the upper construction, anchor isolation, the candidate branch
tree, the exact metric partition, 7,231 noncandidate residual certificates,
negative fail-closed tests, and orbit completeness for all 1,313,024 symmetry
classes. The partition is

\[
1,313,024=1,305,792+7,231+1.
\]

This is a replayable computer-assisted proof claim. It is not a Lean/Coq
formalization and has not yet undergone independent peer review.
