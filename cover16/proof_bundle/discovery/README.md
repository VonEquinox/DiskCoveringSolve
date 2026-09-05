# Optional discovery tools — not an acceptance path

These scripts produced the data checked by the standard-library verification core. They require NumPy, SciPy, mpmath, and numba. Verification itself requires none of these packages. Software versions and numerical optimizer outputs are not trusted by the proof.

Create a **separate** flat research workspace:

```sh
python3 discovery/prepare_workspace.py /path/to/new/research16
cd /path/to/new/research16
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

The supplied candidate and active-system numerical seed are included. `make_root16.py` refines the seed and produces rational midpoint/preconditioner data; `make_local16.py` proposes a rational congruence matrix. Every output must subsequently pass `verify_root16.py` and `anchor_isolation16.py`. The explicit system and active graph are fixed in `system16.py`; changing a candidate graph also requires changing that system, the upper-cover proof and graph interfaces.

To regenerate the combinatorial input tables in this flat workspace, compile `enumerate16.cpp` and run it for `(B,I)=(11,5),(12,4),(13,3),(14,2),(15,1)`. `metric16.py` processes these tables and uses `cycle_metric16.py` to find integer Farkas covering vectors from negative cycles. It writes the metric certificate and ordered residual list. `find_candidate16.py` links the declared candidate face set to that list.

`barrier16.py` is an optional compiled log-barrier solver for a convex per-rod relaxation. `force_fast16.py` uses it before falling back to the SLSQP discovery method in `force16.py`. The solver's success status is not an acceptance criterion: the resulting rationalized forces must satisfy **exact** integer divergence and a strictly positive radius margin.

`build_candidate16.py` writes `candidate16_tree_curvature.json.gz`; this is the discovered tree whose accepted copy is named `core/candidate16_tree.json.gz`. `build_forest16.py` writes one candidate-free residual tree per case under `forest16/`, and `merge_forest16.py` creates the stream `noncandidate16_trees.jsonl.gz`.

All merged or changed data must go through the independent `core/verify_forest16.py` in a correctly staged proof directory and the main verifier. Discovery scripts may fail, stop at search limits or find different trees. None of those events establishes a theorem. The delivered proof uses only the fixed successful certificate files in `core/` and `enumeration/` and their full exact replay.
