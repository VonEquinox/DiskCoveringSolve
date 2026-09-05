# Cover19: nineteen equal disks

This directory contains two related computer-assisted global optimality
proof packages for

\[
r_{19}=1/\sqrt{13},\qquad R_{19}=\sqrt{13}.
\]

Here `r19` is the minimum common radius covering the unit disk, and `R19`
is the corresponding maximum disk radius covered by nineteen unit disks.
The exact construction uses rational pairs representing coordinates in
`Q(sqrt(3))`; there is no numerical root-isolation step.

Both full chains passed fresh replay twice on September 5, 2026. This is a
global lower-bound claim, not a report of numerical local optimization.
Independent external mathematical and software review remains necessary;
the packages are not proof-assistant formalizations.

## Proofs and methods

- [Primary Chinese proof](docs/r19_proof_zh.md): explicit triangular-lattice
  cover, geometric reduction, complete root-face peeling with safe pruning,
  exact force certificates, and global Dirichlet-energy convexity for the
  candidate graph.
- [Alternative Chinese proof](docs/r19_alternative_proof_zh.md): the same
  exact cover under different interior labels, a rational 70-coordinate
  isolation certificate, and a complete candidate-domain branch tree.
- [Audit record](docs/AUDIT_STATUS.md): checked claims, both replay results,
  exact cross-package correspondence, and remaining trust boundaries.
- [Optional discovery sources](alternative_bundle/discovery/README.md):
  candidate search and certificate construction. These numerical programs
  are not used for acceptance and were not rerun for this release.

The two chains share mathematical reductions and some software ancestry.
Their agreement is not two independent external reviews. Neither requires
the optimality claims for n = 11 through 18 as a premise.

## Run verification

Requirements: Python 3.10+ and a C++17 compiler available as `g++`/`c++`.
No Boost or third-party Python packages are needed for exact acceptance.
Allow several hundred MB for inputs and temporary working copies; actual
run time depends on the compiler and machine. Each original master took
roughly 1-2 minutes on the local release machine, with additional time for
the full cross-bundle orbit comparison.

From the repository root:

```bash
make cover19-full PYTHON_BIN="$(command -v python3)"
```

Required final marker:

```text
COVER19 REQUESTED REPLAYS AND EXACT LINKAGE PASSED
```

Other commands:

```bash
make cover19-integrity
make cover19-runner-tests
make cover19-primary-full
make cover19-alternative-full
make cover19-linkage
```

The individual full commands also check cross-bundle linkage, but only replay
the requested master. Integrity, linkage, and synthetic runner tests alone
are not substitutes for both full chains. Never run the verifiers with
`-O`, `-OO`, or `PYTHONOPTIMIZE`.

Retain new reports in a directory that does not already exist:

```bash
python3 -S -B cover19/run_verification.py --report-dir /tmp/cover19-fresh-reports
```

## Results and scope

Both chains regenerate seven families `(B,I) = (12,7), ..., (18,1)`.
After proof-backed necessary-condition cuts, 621,582 rooted survivors form
24,127 dihedral orbit classes. These are not the total unpruned triangulations.
All 24,126 noncandidate types are excluded by exact continuous-domain
certificates. The remaining candidate is handled globally in the primary;
the alternative checks all 527 nodes of its candidate tree, including 28
local leaves. Both report zero unresolved leaves.

The upper construction has 48 triangles and 12 circular caps. Six triangles
are exactly degenerate; the proof uses a nonnegative-orientation winding
argument, not the false claim that all 48 have strictly positive area.
The positive-stress lower-bound subgraph has 42 rods and omits the central
disk center. Its rigidity statement does not classify all optimal covers.

## Files and integrity

- `proof_bundle/`: all 23 files from `r19_verified_proof.zip`, unchanged.
- `alternative_bundle/`: all 66 files from `r19_verified_proof(1).zip`, unchanged.
- `INPUT_MANIFEST.json`: SHA-256 digests for all 89 original files.
- `docs/`: both standalone Chinese proofs and the repository audit record.
- `run_verification.py`: immutable-input checks and fresh temporary-copy replay.
- `verify_bundle_linkage.py`: exact geometry, stress, and full survivor-set comparison.
- `verification/`: fresh reports, logs, environment information, and checksums.

The supplied alternative reports are retained only as original archive
content. The runner deletes reports in its temporary copy before replay,
requires all stages to execute, and exports only newly created reports.
In particular, supplied `reports/master_failure_tests.json` is not reported
as a freshly executed master stage. No manuscript PDF/LaTeX files are added.
