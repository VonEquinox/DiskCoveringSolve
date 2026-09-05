# Covering the Unit Disk with 16 Equal Disks

This package claims the globally optimal equal-disk covering radius

```text
0.308219801898617848870103269221814832006106420717596171183693
 < r_16 <
0.308219801898617848870103269221814832006106420717596171183694.
```

The exact value is `sqrt(t_star)`, where `t_star` is coordinate 74 (zero-based)
of the unique root in the rational box for the 140-variable integer quadratic
system in `proof_bundle/core/system16.py` and `proof_bundle/core/root16.json`.
The proof does not assume optimality of any smaller disk count.

## Proof and method

The [Chinese mathematical proof](docs/r16_proof_zh.md) is byte-identical to the
submitted manuscript. The chain uses a rigorous upper cover, a Voronoi
reduction of arbitrary smaller covers, exhaustive disk-triangulation orbits,
exact angle/Farkas cuts, and integer force certificates on complete angle
domains. Local isolation is used only in the certified candidate leaves.

| Stage | Claimed certificate count |
|---|---:|
| Topology orbits | 53,059,205 |
| Direct metric exclusions | 51,464,643 |
| Farkas exclusions | 1,594,386 |
| Noncandidate residual types | 175 |
| Candidate residual types | 1 |
| Candidate tree nodes | 10,347 |
| Candidate force / local leaves | 5,052 / 122 |
| Noncandidate forest nodes / force leaves | 713 / 444 |

These are global optimality proof claims, not only numerical optimization
results. See [the audit record](docs/AUDIT_STATUS.md) for the independently
replayed scope, provenance and remaining trust boundary.

## Replay

Requirements: Python 3.10+, assertions enabled, a little-endian machine, a
C++17 compiler and Boost.Multiprecision headers. The acceptance chain uses no
third-party Python packages or numerical optimizers. On macOS the repository
runner detects Homebrew's Boost include directory.

From the repository root:

```sh
make cover16-integrity
make cover16-full PYTHON_BIN="$(command -v python3)"
```

The final marker is `COVER16 EXACT MASTER REPLAY PASSED`. All five censuses
are regenerated before independent breadth-first orbit auditing. Expect
several GB of memory and temporary storage, including roughly 1.9 GiB of
binary census tables. The default is one concurrent family; this is a full
replay, not a quick audit. Runtime depends on the machine.

To use two workers and preserve fresh reports in a new directory:

```sh
python3 -S -B cover16/run_verification.py --jobs 2 --report-dir /tmp/cover16-fresh-reports
```

The runner copies inputs into a disposable directory, checks archived hashes
before and after replay, and removes temporary census tables when finished.
The 90 original files, including submitted reports, remain unchanged.

## Layout

- `proof_bundle/`: complete original archive payload, including discovery tools.
- `docs/r16_proof_zh.md`: the original Chinese mathematical proof.
- `docs/AUDIT_STATUS.md`: mathematical review and fresh replay record.
- `verification/`: fresh verification evidence, separate from submitted reports.
- `INPUT_MANIFEST.json`: SHA-256 hashes of all 90 archived files.
- `run_verification.py`: portable, temporary-copy full replay entry point.

This is not a proof-assistant formalization or completed external peer review.
