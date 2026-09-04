# Verifier architecture

## Unified entry point

`run_verification_portable.sh` copies `proof_bundle/` into a temporary working
directory, replaces only historical `/mnt/data` path literals in that copy,
resolves the Python executable to an absolute path, and runs `verify_all.sh`.
The source certificates remain unchanged.

## Certificate dependency graph

```text
full KKT root certificate
├── exact upper candidate
├── exact equilibrium multipliers
└── local transfer bounds

triangulation enumeration
├── Brown labeled-count check
├── stored orbit equality
└── metric/Farkas residual set
    ├── 3788 Farkas exclusions
    ├── candidate orbit linkage (9,2,2547)
    ├── wheel orbit linkage (10,1,1003)
    └── 53 ordinary residual types

candidate certificates
├── exact proposal Kron rows
├── 106136 base leaves
├── 304708 refined nodes
├── 153052 external leaves
├── 1953 local leaves
└── 8D local Hessian/Taylor theorem
```

## Quick mode

Quick mode verifies:

- rational trigonometric primitives and reduction constants;
- full labeled enumeration and Brown counts;
- Farkas residual classification;
- candidate and wheel orbit linkage;
- rational Krawczyk root isolation;
- exact recomputation of six local Kron rows;
- local transfer and 9-gap/8D linkage;
- candidate tree structure;
- upper-complex topology and geometry;
- wheel row validity;
- aggregation of stored full-replay leaf results.

It deliberately does not recompute every leaf.

## Full mode

Full mode additionally recomputes:

- all 106136 candidate base leaves;
- all refined candidate leaves;
- all wheel leaves;
- all 13789 strict leaves for 53 ordinary residual types.

## Trusted base

The exact trusted base is primarily Python integer arithmetic and
`fractions.Fraction`. Trigonometric functions are bounded by rational Taylor
polynomials with explicit remainders. Matrix positivity uses exact rational
`LDL^T`. NumPy/SciPy are dependencies of proposal or support routines, while
accepted proof claims are rechecked by exact or rigorously bounded code.
