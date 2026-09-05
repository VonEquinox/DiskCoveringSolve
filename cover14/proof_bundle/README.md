# Fourteen equal disks covering the unit disk — proof package

The theorem proved in `PROOF_zh.md` is

\[
r_{14}=\sqrt{t_*}=0.331732034276234122557815548802409538270942754829\ldots,
\qquad R_{14}=1/r_{14}=3.014481257988179142001870494819761487186097134743\ldots.
\]

Here `t_star` is coordinate 62 (zero-based) of the unique root of the explicitly constructed 116-variable integer quadratic system in the rational box supplied by `core/root14.json`.

## Replay

From this directory:

```sh
python3 -S -B verify_all.py
```

Requirements: Python 3.10+, a C++17 `g++` or `clang++` compiler, and Boost.Multiprecision headers. No numerical Python package and no network access is needed. Do not use `-O`, `-OO`, or `PYTHONOPTIMIZE`.

The program builds its independent integer orbit auditor, rechecks all four enumeration families, reconstructs the metric partition from the original faces, and replays every force/local leaf. It also recomputes the algebraic-root and matrix certificates and the exact upper cover. Existing success reports are not accepted as evidence. Final success is written to `reports/MASTER_VERIFIED.json`; failure removes any old master success report.

The final partition is

```
1,313,024 = 1,305,792 metric exclusions + 7,231 noncandidate cases + 1 candidate.
```

The candidate tree has 913 nodes: 456 splits, 419 strict force leaves, and 38 anchor-isolation leaves. The noncandidate forest has 16,075 nodes: 4,422 splits and 11,653 strict force leaves. There are no unresolved leaves.

## Contents

`PROOF_zh.md` contains the complete written geometric reduction and proof. `core/` contains the exact verifiers, algebraic root, positive-definiteness certificate, and candidate tree. `enumeration/` contains the four complete orbit lists, independent C++ auditor, metric certificates, exact candidate mapping, and all noncandidate trees. `reports/` contains the fresh replay reports and logs. `discovery/` contains optional numerical-search sources, which are not imported by the proof verifier.

The upper proof explicitly handles six additional equality faces with zero KKT weights through certified symmetry identities. It also handles two collapsed ring triangles by a topological-degree covering argument. Neither issue is hidden behind a numerical tolerance.

## Scope

This is a computer-assisted mathematical proof with explicit written geometric and topological lemmas and independently replayable finite certificates. The lower-bound reduction does not assume that the 11-, 12-, or 13-disk optimum has already been proved. It is not a Lean/Coq formalization and does not assert that external peer review has occurred.
