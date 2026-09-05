# Optional, untrusted numerical discovery

These scripts were used to propose the candidate and exact certificates. None is imported by `verify_all.py`. Verification works without their numerical dependencies. Proposals are not accepted until the standard-library accepting code reconstructs and proves every required inequality.

Create a separate flat workspace (destination must not already exist):

```bash
python3 discovery/prepare_workspace.py /absolute/path/to/new-r18-discovery
cd /absolute/path/to/new-r18-discovery
```

Optional discovery dependencies are NumPy, SciPy, mpmath and Numba. No commercial solver is needed. The initializer copies the certified candidate geometry (`active18_numeric.json`), the accepting graph/system code, the complete survivor list and the candidate map into the workspace. It does not modify the proof package.

Starting from the supplied numerical seed:

```bash
python3 make_root18.py
python3 verify_root18.py
python3 make_local18.py
python3 build_candidate18.py
python3 build_forest18.py 6
python3 merge_forest18.py
```

`make_root18.py` proposes `root18.json`; the separate exact root verifier accepts it. `make_local18.py` proposes a rational congruence then calls the exact matrix verifier. Candidate discovery writes `candidate18_tree_curvature.json.gz`; the corresponding accepted-package filename is `candidate18_tree.json.gz`. Forest discovery writes separate files in `forest18/`; merging writes `noncandidate18_trees.jsonl.gz`. New proposals must undergo the full packaged verification and may differ byte-for-byte with numerical-library versions; no generated optimizer status constitutes a proof. Cached case files are untrusted proposals and are independently verified after merging.

`discover18_B12.py` and `search18.py B iterations` provide numerical geometry searches; `analyze18.py` extracts active constraints and multipliers from `candidate18_B12_numeric.json`. Changing a candidate may require rebuilding the explicit system and all linkages; do not blindly replace a system or root in the proof package. The included `system18.py` is the fixed, explicit graph proved in the manuscript, not an arbitrary active-set generator.

The rigorous archive is self-contained without repeating the discovery stage.
