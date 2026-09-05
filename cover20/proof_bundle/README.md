# n=20 equal disks covering a unit disk

The exact proposed theorem and complete mathematical reductions are in `PROOF_zh.md`.
The answer is the radius coordinate of the specified 198-variable integer quadratic system, not a floating-point optimizer output.

## From-source replay

```sh
python3 -S -B verify_all.py
```

Requires Python 3.10+ and a C++17 compiler (`g++` or `clang++`). Set `CXX` to select another compiler. The accepting path uses only the Python standard library and standard C++17. No Boost, NumPy, SciPy, MPMath, NumPy/Numba caches, solver, or network is required. Never use `-O` or `-OO` with Python.

The program regenerates all eight safely pruned exhaustive topology families twice, compares complete canonical survivor sets, checks the packed representatives, replays both rational contraction boxes and the complete upper cover, checks the anchor-only matrix theorem, and verifies every case and every leaf. Existing success reports are never used instead of recomputation.

`COVER20_WORKERS=4` is the default. Any integer from 1 to 16 changes task parallelism only. The decompressed recursion outputs are placed in `_runtime/` and can be removed after a run. Allocate several gigabytes of free disk space and memory; the raw proof forest is streamed, not loaded all at once.

The compressed topology representatives and force trees are proof inputs. Do not omit them. They are checked against freshly generated data and exact inequalities. The force format stores cycle-edge forces; tree-edge forces are reconstructed exactly and all node balances are recomputed independently.

After success, `reports/MASTER_VERIFIED.json` records all input hashes and exact counts. Any failure deletes that master success file. Individual logs and reports are retained for diagnosis.

`discovery/`, when present, contains optional numerical proposal scripts. It is not imported by the verifier. Numerical discovery uses additional dependencies and is not part of the proof's trust base.

## Scope

This is a computer-assisted proof with written convex-geometric/topological reductions. It is not a Lean/Coq formalization and has not been claimed to have undergone third-party peer review. The theorem concerns the optimal radius, not a classification of all optimal configurations. The known numerical scale R=3.692+ was previously recorded by Erich Friedman as a Jeremy Tan (2018) construction.

The distribution ZIP intentionally excludes old success reports and `_runtime/`. A fresh run creates all reports itself. Published replay reports are supplied as separate companion downloads. No success flag is a proof input.
