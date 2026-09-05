# Covering the Unit Disk with 17 Equal Disks

The package claims the globally optimal equal-disk covering radius

```text
0.298531589256428059869843774423980703734387603913740470124900
 < r_17 <
0.298531589256428059869843774423980703734387603913740470124901.
```

The exact value is `sqrt(t_star)`, where `t_star` is coordinate 82 (zero-based)
of the unique root in the specified rational box of the 157-variable integer
quadratic system. See `proof_bundle/core/system17.py` and `proof_bundle/core/root17.json`.
The proof does not assume optimality for any smaller number of disks.

## Proof and method

The [Chinese proof](docs/r17_proof_zh.md) gives the geometric reduction,
exhaustive root-face peeling, safe partial-state pruning, upper cover,
integer force inequalities, local isolation and complete domain partition.
It is proof documentation, not an arXiv submission package.

Unlike Cover16, this package does not first list every unpruned completed
triangulation. It exhausts root-face choices and rejects partial states only
using necessary four- and six-rod angle inequalities. The surviving completed
graphs are independently regenerated using different region/root/label
orders, BFS canonicalization and full Floyd-Warshall closure. Their complete
sets are compared, not merely their counts.

| Family | Primary partial states | Independent partial states | Surviving orbits |
|---|---:|---:|---:|
| B11 I6 | 161,886 | 102,179 | 489 |
| B12 I5 | 177,422 | 242,611 | 631 |
| B13 I4 | 90,653 | 360,159 | 200 |
| B14 I3 | 44,600 | 477,303 | 0 |
| B15 I2 | 26,820 | 429,475 | 0 |
| B16 I1 | 15,293 | 204,649 | 24 |
| Total | 516,674 | 1,816,376 | 1,344 |

**1,344 is the number of surviving orbits, not the total number of all
unpruned triangulations.** The written peeling induction and safe-pruning
argument are essential to completeness.

All 1,343 noncandidate survivors are excluded by 2,709 force leaves. The
candidate's entire angle domain has 1,357 force leaves and 37 local leaves;
there are zero unresolved leaves. The upper cover uses 44 triangles and
12 circular caps. The local theorem applies to 63 active rods at anchor
radius 27/250. See [the audit record](docs/AUDIT_STATUS.md) for checked scope
and trust limits.

## Replay

Python 3.10+ with assertions enabled and a C++17 compiler are required.
No Boost, third-party Python packages, numerical optimizer or network is
required for acceptance.

```sh
make cover17-integrity
make cover17-full PYTHON_BIN="$(command -v python3)"
```

The final marker is `COVER17 EXACT MASTER REPLAY PASSED`. Every run rebuilds
both C++ programs and regenerates all six surviving sets. The repository
runner removes stale success reports in its temporary copy and leaves all
46 archived files untouched. Temporary compilation/search files are deleted
afterward.

To retain fresh stage reports in a new directory:

```sh
python3 -S -B cover17/run_verification.py --report-dir /tmp/cover17-fresh-reports
```

The master includes 22 unpruned small-count regressions and 30 fail-closed
tests. The additional repository test extends unpruned counting to four,
five and six interior vertices:

```sh
make cover17-extra-tests PYTHON_BIN="$(command -v python3)"
```

## Contents and limits

- `proof_bundle/`: all 46 original submitted files, including supplied reports.
- `docs/r17_proof_zh.md`: byte-identical standalone Chinese proof.
- `docs/AUDIT_STATUS.md`: mathematical review, provenance and fresh-run record.
- `verification/`: fresh replay reports, kept apart from submitted reports.
- `INPUT_MANIFEST.json`: all archived file hashes.
- `run_verification.py`: temporary-copy complete replay entry point.
- `test_peeling_counts.py`: additional unpruned counting regressions.

No PDF, LaTeX manuscript or arXiv submission files are included. Numerical
discovery programs were not supplied in this archive; they are not accepting
dependencies. The complete candidate, finite certificates and verification
sources are included. External peer review and proof-assistant formalization
are not claimed.
