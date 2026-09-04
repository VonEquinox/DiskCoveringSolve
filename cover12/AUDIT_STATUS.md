# Cover12 audit status

Audit date: September 4, 2026.

The bundle passed a fresh local execution of `cover12_master_verify.py` using
Python 3.13.5, NumPy 2.5.2, SciPy 1.18.0, and SymPy 1.14.0. The replay checked
all eight master modules and returned code 0 in approximately 30 seconds.

The following links were checked during the audit:

- the fourteen-variable rational Krawczyk isolation certificate;
- symbolic contact, equilibrium, normalization, and anchor-stationarity identities;
- exact threshold and trigonometric comparisons;
- the 31-face upper-covering complex;
- validity, canonicality, orbit sizes, and Brown totals for all six triangulation families;
- exact Farkas certificates for 38,076 excluded orbits;
- exact set equality between the 132 residual orbits, the 131 noncandidate energy certificates, and the single candidate orbit;
- exact Kron reduction and every stored noncandidate branch leaf;
- positive rational `LDL^T` pivots for the candidate global-convexity certificate.

The mathematical exposition was also strengthened in two places before this
repository release: the boundary no-return argument now deletes precisely the
cells owning open subarcs of the intervening gap, and the passage from the
restricted Voronoi subdivision to a simple disk triangulation is stated and
proved as a separate proposition.

No contradiction, failed exact check, or certificate-linkage mismatch was
found. This audit is not a substitute for independent peer review or a formal
proof assistant verification. The appropriate status is a reproducible
computer-assisted proof claim with a passing exact verifier.
