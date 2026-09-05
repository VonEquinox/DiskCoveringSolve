# Optional, untrusted numerical proposal generation

Nothing in this directory is imported by `verify_all.py`. This directory is supplied so that the numerical force and matrix proposals can be rediscovered; **neither optimizer status nor a numerical objective is an accepting condition**.

The discovery programs need NumPy, SciPy, and Numba. They are not needed for the proof replay. Minor optimizer/platform differences can produce different, equally valid trees; bit-for-bit reproduction of compressed proposals is not asserted.

Create a separate workspace (the path must be new or empty and outside the proof package):

```bash
python3 discovery/prepare_workspace.py ../r19_discovery
cd ../r19_discovery
python3 make_local19.py
python3 build_candidate19.py
python3 build_forest19.py 4
python3 merge_forest19.py
```

`make_local19.py` proposes a rational congruence using a floating Cholesky factor and immediately submits it to the exact matrix verifier. It does not change the exact weights, geometry, or theorem constants.

`build_candidate19.py` uses the logarithmic-barrier/SLSQP routines to propose forces and split decisions. Its output is `candidate19_tree_curvature.json.gz`. `build_forest19.py` proposes per-case noncandidate trees under `forest19/`. `merge_forest19.py` checks that every expected noncandidate case is present and constructs `noncandidate19_trees.jsonl.gz`.

`find_candidate19.py` is an optional graph-isomorphism discovery tool. The proof does not trust its search; it directly checks the saved permutation.

The topology survivor manifest is an input to these proposal tools, but is **not** an unverified input to the proof: the main verifier recompiles both C++ recursions, regenerates all survivors, and compares every graph with that manifest. A replacement proposal must be checked in a separate copy of the proof package with the complete `verify_all.py`, with all corresponding exact interfaces preserved. A new valid tree may have different counts, requiring an independently justified update of the descriptive count assertions; the exact leaf and completeness tests must never be relaxed.
