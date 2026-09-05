# Exact verification bundles for fifteen disks

The included computer-assisted proof claims

\[
r_{15}=\sqrt{t_*}=0.318142930859262836359491181649935005186663866157\ldots.
\]

Here `t_star` is coordinate 76 (zero-based) of the unique root in the rational
box for the 148-variable integer quadratic system. Both bundles define this
system explicitly in `core/system15.py` and supply `core/root15.json`. Their
labels differ; `verify_bundle_linkage.py` identifies the roots by exact graph
relabeling and rotation into the primary uniqueness box.

The strict certified interval is

```text
0.318142930859262836359491181649935005186663866157651845538266
 < r_15 <
0.318142930859262836359491181649935005186663866157651845538267.
```

## Contents

- `proof_bundle/`: contents of `r15_verified_proof.zip`;
- `alternative_bundle/`: contents of `r15_verified_proof(1).zip`;
- `docs/r15_proof_zh.md`: separately submitted Chinese manuscript, identical
  to `alternative_bundle/PROOF_zh.md`;
- `docs/AUDIT_STATUS.md`: review scope, provenance, and comparison;
- `verification/`: fresh replay logs, master reports, and exact root linkage;
- `INPUT_MANIFEST.json`: all 146 archived files, including submitted reports
  and the standalone manuscript;
- `run_verification.py`: archive verification and temporary-copy full replay.

## Replay

Requirements: Python 3.10+, a C++17 compiler, Boost.Multiprecision headers, and
a little-endian machine for the alternative binary format. Allow several GB
of free disk space and memory, and several minutes for regeneration and audit.
No numerical Python library is required for verification.

On macOS, install Boost with `brew install boost`. The runner detects its
Homebrew include directory. On Debian/Ubuntu, install `g++ libboost-dev`.

From the repository root:

```bash
make cover15-full PYTHON_BIN="$(command -v python3)"
```

This runs both complete verifier chains, then independently rechecks both root
certificates and their exact correspondence. Its final marker is:

```text
COVER15 REQUESTED REPLAYS AND ROOT LINKAGE PASSED
```

Individual replays and archive checks:

```bash
make cover15-primary-full
make cover15-alternative-full
make cover15-integrity
```

The runner uses temporary copies and regenerates the entire census. Original
archive reports and hashes stay fixed. Numerical discovery scripts are
included for research, but are not used by the accepting verifiers.

## Scope

Each chain audits 11,950,884 topology orbits and checks the partition
`11,907,870 + 43,013 + 1`. Both use a 76-dimensional anchor-isolation argument
and an exact upper cover with 26 vertices, 64 edges and 39 triangles.

These packages share some source and mathematical arguments. Successful
replays are evidence for their finite certificate arithmetic, not two
independent external peer reviews or a proof-assistant formalization.
