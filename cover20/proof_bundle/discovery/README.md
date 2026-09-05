# Optional numerical discovery, outside the accepting path

These are proposal tools, not proof verifiers. Their residuals, success flags,
floating-point objective values, and numerical Hessians are not evidence of a
theorem. The complete accepting path is `../verify_all.py` and does not import
this directory.

To work in an empty separate directory:

```sh
python bootstrap.py /path/to/new-empty-workspace
cd /path/to/new-empty-workspace
```

The numerical tools use NumPy, SciPy, MPMath and Numba. `bootstrap.py` copies
accepted data and proposal kernels into that workspace, and expands the packed
masks into NumPy arrays. It refuses to write into the proof directory.

The archived `active20_numeric.json` provides a seed for `make_root20.py`.
`make_root20.py` proposes a high-precision root and rational inverse;
`verify_root20.py` must independently accept the new root.
`make_local20.py` proposes a triangular congruence and then calls its exact
acceptor. `verify_upper20.py` checks the entire-disk upper cover.

`search20.py 13 120 1` runs numerical candidate exploration. It is not exhaustive
and it may converge to a different local solution. `analyze20.py` builds a
numerical active-set seed for the archived topology; a different active set
requires rebuilding the explicit polynomial system rather than silently using
the old one.

`build20.py candidate` proposes the candidate tree.
`build20.py all 4` proposes noncandidate trees in atomic chunks of 500 cases.
Every stored proposal leaf is checked exactly during discovery, but an
independent full accepting replay is still required. Any failed chunk is
reported separately and is not saved as a completed chunk. A successful
numerical run does not authorize a claim of global completeness.

The prepared mask arrays already come from the full accepted topology inputs.
To verify that those inputs are exhaustive, use the from-source proof replay,
not this optional discovery workspace.
