# Cover13 audit status

**Independent replay date:** September 5, 2026 (Asia/Shanghai)

## Result

A fresh replay was executed from an independent extraction of the submitted
ZIP after verifying its internal ZIP integrity and all 28 SHA-256 manifest
entries.  The replay returned code 0 and ended with all stages marked `PASS`.

The generated aggregate report states:

- `verified: true`;
- `finite_certificates_verified: true`;
- `enumerated_cases_complete: true`;
- 308,198 total topology orbits;
- 306,437 metric exclusions;
- 1,760 noncandidate residual exclusions;
- one candidate orbit;
- 1,831 candidate tree nodes with 837 force leaves and 79 local leaves;
- 2,908 noncandidate tree nodes with 2,334 force leaves.

The independent replay log is
`../verification/independent_full_replay.log`, and its aggregate JSON is
`../verification/independent_MASTER_VERIFIED.json`.

## Environment note

On this macOS/Homebrew system, Boost is keg-only.  The original master verifier
compiled successfully after exposing `$(brew --prefix boost)/include` through
`CPLUS_INCLUDE_PATH`.  The repository wrapper performs this detection
automatically.

## Scope and limitations

The replay establishes that the archived finite exact-arithmetic checks pass
and that the enumeration audit agrees with the written recurrence counts.  It
does not substitute for independent review of the analytic, geometric, and
PL-topological lemmas in `cover13_proof_zh.md`, nor does it constitute a
Lean/Coq formalization or peer review.
