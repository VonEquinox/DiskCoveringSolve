# Cover20: twenty equal disks

This package claims the globally optimal common radius for covering the unit
disk with twenty equal closed disks:

\[
r_{20}=0.27084817843997736519402387625893583587132712\ldots.
\]

The exact value is the square root of the radius-squared coordinate of the
specified 198-variable integer quadratic system in its certified rational
root box, not the displayed decimal or an optimizer's output.

Both the direct ZIP replay and the repository reconstruction/full replay
passed locally on September 5, 2026, with zero unresolved leaves.

- [Chinese proof from the ZIP](docs/r20_proof_zh.md)
- [Audit record and proof scope](docs/AUDIT_STATUS.md)
- [Optional discovery methods and source](proof_bundle/discovery/README.md)
- [Fresh verification evidence](verification/README.md)

Only the proof supplied inside `r20_verified_proof.zip` is used. The separate
Markdown submission describes different root data, local constants and force
formats; it is not included as a second verified proof.

## Full replay

From the repository root:

```bash
make cover20-full PYTHON_BIN="$(command -v python3)"
```

Requirements: Python 3.10+ and a C++17 compiler (`g++` or `clang++`; `CXX` may
select another executable). The accepting path needs no third-party Python
packages, Boost, network, or numerical optimizer. Allow several GB of free
disk space and memory. The direct release audit took about 7.5 minutes on the
recorded local machine; other machines may take substantially longer.

Use a lower worker count to reduce parallel resource use, and a new directory
to retain fresh reports:

```bash
python3 -S -B cover20/run_verification.py --workers 2 --report-dir /tmp/cover20-fresh
```

The default is four workers. Worker count changes task parallelism only, not
the accepting inequalities, search completeness, or input data. Every full
run regenerates both exhaustive searches for all eight families and checks
every continuous case. No cached-success or partial-enumeration mode is used.

Required final marker:

```text
COVER20 EXACT MASTER REPLAY AND INPUT RECONSTRUCTION PASSED
```

Lightweight checks, not substitutes for full replay:

```bash
make cover20-integrity
make cover20-runner-tests
```

## Lossless large-file storage

The original `enumeration/noncandidate20_trees.jsonl.gz` is 197,440,381 bytes,
too large for an ordinary GitHub Git blob. It is stored as three consecutive
byte ranges in `proof_parts/`, without recompression or alteration.

`INPUT_MANIFEST.json` records the original SHA-256 of every one of the 51 ZIP
files. `STORAGE_MANIFEST.json` records the ordered part names, sizes and hashes.
Integrity mode hashes the concatenated byte stream and checks the original
digest. Full mode reconstructs the original file in a temporary copy, verifies
all 51 files, then invokes the unchanged original `verify_all.py`.

No original ZIP file is discarded. Fifty files are stored directly under
`proof_bundle/`; the remaining logical file is restored from the parts.
No LFS setup or external download is needed. Do not run `proof_bundle/verify_all.py`
directly from the stored tree: use the repository runner so the missing large
logical file is restored first. The original upstream README describes a
fully reconstructed bundle, not this split storage layout.

## What is checked

The upper bound uses 51 positively oriented triangles and 13 circular caps.
For the lower bound, ordinary Voronoi reduction and exhaustive root-face
peeling cover eight families `(B,I) = (12,8), ..., (19,1)`. Necessary-condition
pruning leaves 562,169 orbit classes; this is not the unpruned topology count.
All 562,168 noncandidate cases are checked through exact integer force
certificates. The candidate's 6,171-node complete tree combines force
exclusions with a 102-coordinate anchor-isolation theorem of radius `11/100`.
The full chain also runs 33 unpruned count regressions and 43 rejection tests.

An additional repository checker independently reconstructs all 198 residuals
and 39,204 Jacobian entries from squared distances and checks the ZIP proof's
face tables, root hash and displayed centers against the actual system.

This is a computer-assisted global optimality proof claim, not just a local
search result. Written geometry, topology, analysis and enumeration induction
remain part of the trust boundary. Independent external review is still
needed; the package is not a Lean/Coq formalization and does not classify all
optimal covers. No arXiv-style PDF/LaTeX manuscripts are added.
