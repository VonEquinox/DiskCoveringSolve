# Exact certificate package for the 18-disk covering theorem

中文完整证明见 `PROOF_zh.md`。

## Fresh verification

```bash
cd r18_verified_proof
python3 -S -B verify_all.py
```

Requires Python 3.10+ and a C++17 compiler (`g++` or `clang++`) on PATH. No third-party Python packages, Boost, numerical optimizer, or network are required for verification. Do not use `-O`/`-OO` or set `PYTHONOPTIMIZE`. Run from a writable extracted directory; allow space for newly generated topology outputs in `_runtime/`.

The command regenerates both independent exhaustive, safely pruned topology searches for all seven families. It then verifies equality of complete survivor sets, the exact algebraic root, the actual covering, the quantitative anchor isolation matrix, all angle-tree leaves, 29 unpruned counting tests, and 32 fail-closed rejection tests. A stored `verified` flag never substitutes for a fresh check.

Success is written to `reports/MASTER_VERIFIED.json`; details are in `reports/FULL_REPLAY.log` and per-stage logs. Any failure deletes the master success report. Input hashes must remain unchanged during the run.

## Scope

The exact optimum is the radius coordinate of the isolated 174-variable integer quadratic KKT root. Display approximations are r = 0.2901677176405667291131545..., R = 3.4462827503048036561495422....

There are 16,220 surviving topology orbits AFTER exhaustive safe partial-graph cuts, not 16,220 total unpruned triangulation orbits. The proof covers the pruned extensions by a written safety/completeness argument, and all survivors by 16,219 complete exclusion trees and one complete candidate tree. The latter uses an anchor-only quantitative theorem on a radius 11/100 neighborhood. No smaller-n optimality theorem is assumed.

This is an exact-arithmetic computer-assisted proof with written geometric and finite-recursion reductions. It is not a Lean/Coq formalization and is not a claim of completed external peer review.

## Files

`core/` holds accepting code and algebraic/matrix/candidate-tree certificates. `enumeration/` holds the two recursive C++ programs, the complete surviving face lists, candidate mapping, and all other angle trees. `_runtime/` is regenerated scratch output and need not be distributed. `discovery/` is optional numerical discovery code and is NOT trusted by verification. `SHA256SUMS` provides distribution checksums.

All supplied Python accepts arbitrary working directories because verification inputs are located relative to the package. Discovery programs are separate; see their own README for preparing a flat numerical workspace.
