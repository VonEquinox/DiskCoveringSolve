# Exact computer-assisted proof for 19 congruent disks

The theorem in `PROOF_zh.md` is

    r_19 = 1/sqrt(13),   R_19 = sqrt(13).

Run from this directory:

```bash
python3 -S -B verify_all.py
```

Requirements: Python 3.10 or later (standard library only), and a C++17 compiler (`g++` or `clang++`). No Boost, numerical solver, internet connection, or third-party Python package is needed. Assertions must remain enabled; `-O`, `-OO` and `PYTHONOPTIMIZE` are rejected.

The replay rebuilds all exact geometric data from a 19-point lattice, certifies a 70-by-70 rational matrix, recompiles and reruns both independent root-face recursions on B=12,...,18, compares every surviving graph, and replays every force/local leaf. It also reruns 29 unpruned counting tests and 36 rejection tests. It does not trust old success reports or use numerical optimizer status as proof.

The final acceptance report is `reports/MASTER_VERIFIED.json`; a failed run removes any old report at that path. Exact certificates occupy roughly 100 MB compressed. The full replay takes several minutes depending on hardware; intermediate JSON and compilation outputs are written into `_runtime/`.

## Important mathematical distinctions

The 24,127 graph orbits are those remaining **after** a proved exhaustive recursion with safe partial-graph cuts. They are not the count of all unpruned triangulations. Every discarded extension is covered by the necessary-condition pruning theorem in the manuscript.

The 42 positive-weight rods omit the central disk center. The local theorem fixes only the 36-node stressed graph, not all nodes of the full 55-node auxiliary graph. This is sufficient to exclude a strictly smaller radius and is explicitly reflected in the graph embedding verifier.

No Krawczyk root certificate is needed in this case: all candidate coordinates are explicit in Q(sqrt(3)), and their squared distances and stationarity equalities are rational after a documented coordinate scaling. Six covering triangles are exactly degenerate and are treated as such.

This is a computer-assisted proof with written mathematical reductions. The package does not claim external peer review or Lean/Coq formalization. The trust boundary is specified in the manuscript.

`discovery/` contains optional, untrusted numerical certificate-proposal tools. They are not imported by any accepting routine. Their dependencies are unnecessary for verification.
