# Nineteen congruent disks covering a disk

The enclosed computer-assisted theorem is

\[
r_{19}=1/\sqrt{13},\qquad R_{19}=\sqrt{13}.
\]

`PROOF_zh.md` gives the written geometric reduction, exact lattice construction,
42-rod positive stress, globally convex boundary energy, exhaustive rooted-face
recursion, and exact force-certificate lemma. The input certificates are not
accepted by reading their old success flags.

## Full replay

From this directory, run:

```sh
python3 -S -B verify_all.py
```

Requirements: Python 3.10 or later and either `g++` or `clang++` supporting C++17.
The acceptance path uses only the Python standard library and the C++ standard
library. No NumPy, SciPy, linear programming solver, Boost, plantri or network is
needed. Do not use `-O`, `-OO` or `PYTHONOPTIMIZE` for Python.

Every invocation recomputes the exact candidate construction and lower matrix,
regenerates all seven pruned combinatorial families, independently re-enumerates
and compares the exact survivor sets, replays every force leaf, checks the full
candidate/noncandidate partition, and runs rejection regressions. There is no
partial or skip mode.

A successful replay writes `reports/MASTER_VERIFIED.json`. Old master reports are
removed before verification and on failure. Proof inputs are hashed before and
after the run. Stage logs are written to `reports/` and stage certificates are
recomputed in `core/`.

## Main inputs

- `core/geometry19.py`: explicit rational coordinates in the representation
  `(x,y) = (x,sqrt(3)*y)` and all face, edge and weight definitions.
- `core/verify_candidate19.py`: continuous covering upper bound and global
  candidate lower bound, both recomputed with exact rational arithmetic.
- `enumeration/peel_enum.cpp` and `peel_audit.cpp`: independent complete rooted
  face recursions, using different region orders, labels, roots, canonical
  traversals and difference-constraint algorithms.
- `enumeration/B*_I*.txt`: finite survivor manifests, fully regenerated.
- `enumeration/metric19_residuals.json`: exact union of those manifests.
- `enumeration/candidate19_map.json`: explicit candidate graph identification.
- `enumeration/noncandidate19_trees.jsonl.gz`: all 24,126 remaining angle trees.
- `core/exact_force19.py` and `verify_forest19.py`: integer force balance,
  outward arc support, norm bounds, tree coverage and set-completeness checks.

The full 19-center candidate has zero stresses on many active equal-length
constraints. The lower proof deliberately uses only its 42 positive rods and
eliminates their free vertices exactly. It does not claim that all 19 centers
are isolated by a positive-definite full-coordinate KKT system. There are no
candidate local leaves and no numerical root-isolation oracle.

## Scope

This is a computer-assisted proof with written mathematical reductions. Its two
independent in-package enumerations are not an external peer review or a
Lean/Coq formalization. The theorem concerns the optimal radius; it does not
classify all optimal coverings.

The public construction value `sqrt(13)` was already recorded in Erich
Friedman's *Circles Covering Circles*, attributed there to Jeremy Tan (2018).
This package does not claim first discovery of that construction.
