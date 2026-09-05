# Exact computer-assisted proof: 17 congruent disks covering the unit disk

The exact optimum is the square root of the designated coordinate of the unique
root of `core/system17.py` in `core/root17.json`. Read `PROOF_zh.md` for the full
geometric, topological, and exhaustive-peeling reductions.

```
0.298531589256428059869843774423980703734387603913740470124900
  < r_17 <
0.298531589256428059869843774423980703734387603913740470124901
```

## Fresh verification

From the extracted directory:

```sh
python3 -S -B verify_all.py
```

Requirements: Python >=3.10 and `g++` or `clang++` with C++17 support. The accepting
pipeline uses only the Python standard library and two bundled C++ source files.
No network, numerical optimizer, Boost, plantri, or external Python package is
required. Do not use `-O`, `-OO`, or `PYTHONOPTIMIZE` to disable Python assertions.

The default command recompiles both enumeration implementations in a temporary
directory, regenerates all six families, and compares the complete surviving
sets. It then replays all exact force/local leaves. There is no flag that skips
an accepting mathematical stage and still writes the master success report.

The primary search uses stack-based root-face peeling, ascending fresh labels,
DFS canonicalization, and incremental integer difference-bound closure. The
independent audit uses FIFO regions, different root edges, descending fresh
labels, BFS canonicalization, and full integer Floyd-Warshall reconstruction.
Necessary four/six-rod inequalities safely prune incomplete triangulations.
The proof of exhaustive coverage of all potential counterexamples is in the
written proof; agreement of counts alone is not the completeness premise.

The successful result is `reports/MASTER_VERIFIED.json`. The program removes an
old master report before starting and on any failure. Each stage's actual
stdout/stderr is saved in `reports/`. Standalone module summaries are not a
substitute for running the full master command. Proof inputs are hashed before
and after execution and all stage interfaces are checked.

## Expected exact finite result

There are 1,344 surviving topology orbits: 489, 631, 200, 0, 0, 24 for B=11..16.
All 1,343 noncandidate orbits are excluded by 2,709 exact force leaves. The unique
candidate tree has 1,357 exact force leaves and 37 anchor-isolation leaves.
There are zero unresolved leaves. All other branches of the exhaustive
combinatorial search are rejected by necessary inequalities at partial states.
These are not counts of every unpruned completed triangulation.

The candidate root has 157 variables, 63 strictly positive rod weights, and 13
active center faces. The upper cover is verified with 44 triangles and 12 caps.
The exact anchor-only isolation radius is 27/250. The supplementary regression
suite uses 22 unpruned small-case counts and 30 fail-closed corruption tests.

## Attribution and scope

The public covering table by Erich Friedman attributes the approximate R=3.349+
17-disk configuration to Jeremy Tan (2018). This package does not claim first
discovery of that numerical candidate. Its mathematical claim is global
optimality of the precisely specified algebraic radius.

Discovery optimization is not an accepting dependency; the complete exact
candidate, rational inverse, matrix congruence, and force certificates are
included. All finite verification inputs are bundled.

This is not Lean/Coq formalization or a claim of third-party peer review. The
written reductions, Python exact arithmetic, and the C++ compiler remain in the
trust boundary. No classification of all equal-radius optimizers is asserted.
