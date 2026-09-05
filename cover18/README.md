# Cover18: covering a disk with eighteen equal disks

Two related computer-assisted global-optimality proof packages for
`r_18 = 0.29016771764056672911315452560106171830628088...`.
The exact value is the square root of the radius coordinate of a specified
isolated 174-variable integer polynomial root, not a rounded numerical fit.
Both packages are self-contained with respect to Covers 11-17.

- [Primary Chinese proof](docs/r18_proof_zh.md)
- [Alternative Chinese proof](docs/r18_alternative_proof_zh.md)
- [Review, exact linkage and replay evidence](docs/AUDIT_STATUS.md)
- [Optional discovery methods](alternative_bundle/discovery/README.md)

These are global proof claims with written geometric and recursive
completeness arguments, not merely locally optimized configurations. They
remain subject to external mathematical and software review and are not
Lean/Coq formalizations. No manuscript PDF or LaTeX is included.

## Replay

From the repository root, with Python 3.10+ and a C++17 compiler:

```sh
make cover18-integrity
make cover18-runner-tests
make cover18-full PYTHON_BIN="$(command -v python3)"
```

No Boost or third-party Python packages are needed. Both variants recompile
their enumerators, regenerate all seven families, check all surviving graphs
and continuous branches, and execute their corruption/count regressions.
The repository runner additionally rechecks both roots and identifies their
exact values by rational interval transport, then compares all survivor sets.

The final marker is:

```text
COVER18 REQUESTED REPLAYS AND EXACT LINKAGE PASSED
```

Individual chains: `make cover18-primary-full` and
`make cover18-alternative-full`. These also run the exact linkage check.
To keep fresh reports:

```sh
python3 -S -B cover18/run_verification.py --report-dir /tmp/cover18-new-reports
```

Choose a directory that does not exist. The runner never regenerates outputs
inside the original payloads. The alternative upstream entry point can leave
an old report when it rejects `-O`; the repository runner removes all old
reports in a disposable copy before launching normal, assertion-enabled
verification. Check the current exit status, not a stored `verified` flag.

## Scope

| Check | Primary | Alternative |
| --- | ---: | ---: |
| Geometric variables / KKT variables | 90 / 174 | 90 / 174 |
| Active rods / center faces | 72 / 22 | 72 / 22 |
| Upper-cover triangles / circular caps | 46 / 12 | 46 / 12 |
| Surviving topology orbits | 16,220 | 16,220 |
| Noncandidate cases | 16,219 | 16,219 |
| Candidate tree nodes | 1,219 | 2,323 |
| Noncandidate forest nodes | 23,531 | 30,505 |
| Local anchor radius | 7/60 | 11/100 |
| Unpruned count regression families | 26 | 29 |
| Corruption rejection tests | 32 | 32 |
| Unresolved leaves | 0 | 0 |

The orbit count is the complete survivor set after safe partial-graph cuts,
not the number of all unpruned triangulations. Written no-return, duality,
padding and root-face recursion arguments are essential parts of the proof.

## Files

`proof_bundle/` preserves all 31 original files from `r18_verified_proof.zip`.
`alternative_bundle/` preserves all 76 files from `r18_verified_proof(1).zip`,
including its supplied reports and optional discovery scripts.
`INPUT_MANIFEST.json` hashes all 107 original files; `docs/` holds the two
unchanged proofs and this repository's audit. New replay evidence is kept
separately under `verification/`. Do not combine individual certificates
from the two variants.

The optional discovery tools require additional numerical packages documented
in their README. Their execution is not required for acceptance, and this
release does not claim a fresh rerun of the numerical discovery process.
