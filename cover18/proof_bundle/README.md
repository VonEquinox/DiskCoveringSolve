# Equal-disk covering of the unit disk: n=18

This package contains the written Chinese proof and a fail-closed exact-arithmetic verifier for the radius defined by `core/system18.py` and `core/root18.json`. The optimal-radius claim is for arbitrary configurations, not just the candidate topology. It does not assume optimality for smaller n. It is not a Lean/Coq formalization or a claim of external peer review.

## Replay

```bash
python3 -S -B verify_all.py
```

Requirements: Python 3.10+ and a C++17 compiler (`g++` or `clang++`). Only the Python and C++ standard libraries are used. No external Python packages, Boost, numerical optimizer, or network access are needed. Do not disable Python assertions (`-O`, `-OO`, `PYTHONOPTIMIZE`).

Every default run recompiles both enumeration programs, regenerates all seven surviving topology families, independently audits the exact sets, executes 26 unpruned count regressions, rebuilds the 174-variable root system, checks the continuous upper cover and anchor-only isolation matrix, replays every angle-tree leaf, and performs 32 rejection tests. The final report is `reports/MASTER_VERIFIED.json`; stage logs are under `reports/`. Old master reports are removed at start and upon failure. Input hashes are compared before and after the run.

## Contents

- `PROOF_zh.md`: full mathematical reduction, lemmas, system definition, and closure argument.
- `core/`: explicit KKT system, rational root, interval/trigonometric and matrix checks, force certificates, branch replay, and negative tests.
- `enumeration/`: two structurally different exhaustive root-face recursions, seven complete surviving manifests, canonical residual records, candidate map, and the noncandidate certificate forest.
- `verify_all.py`: the accepting orchestrator.

Discovery used floating-point optimization and linear programming only to propose certificates. Neither solver outputs nor previously stored success flags are trusted during acceptance. The exact witness is the polynomial system plus its rational isolation box; printed decimals are certified display intervals. The proof establishes the minimum radius, not a classification of all optimal configurations.
