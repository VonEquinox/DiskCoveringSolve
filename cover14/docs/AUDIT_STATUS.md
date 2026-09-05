# Cover14 audit status

## Source package

- Input archive: `r14_verified_proof.zip`
- Archive SHA-256:
  `e53f62b1bad4e6e601145ccc1f99d283e1552cbe29117a26978a6f67cb6294c0`
- Archive integrity: `unzip -t` passed.
- Original `SHA256SUMS.txt`: all 61 listed files passed before replay.
- Original submitted master SHA-256:
  `08a087c9014e25e2b557a3f9d7a287024458490f71f2722999cd4d507dd6ce0a`

The original package reports are preserved under
`../verification/submitted_package_reports/`. The repository manifest covers
the 37 fixed non-report files and excludes regenerated reports and timing data.

## Independent replay

On September 5, 2026, the package was freshly extracted and run with:

```bash
python3 -S -B verify_all.py
```

On macOS the Homebrew Boost include directory was supplied through
`CPLUS_INCLUDE_PATH`. Existing success reports were overwritten by the master
runner. The replay returned code 0 and ended with:

```text
R14 COMPLETE CERTIFICATE CHAIN VERIFIED
```

The independent outputs are stored as:

- `../verification/independent_MASTER_VERIFIED.json`
- `../verification/independent_full_replay.log`

Their SHA-256 digests are, respectively:

```text
f8d3fcf4a1e3d4dc6bcbaed762faca84581e32466d842e7ad220695f5370c2d2
8cbe409b76cd84047fd1c5726bca20202f4eebc0aeac768c4ae6546a12e7b6c4
```

## Replayed claims

- the 116-variable root isolation and strict decimal interval;
- the exact upper cover, including six zero-weight equality faces and two
  collapsed ring triangles;
- the radius-`1/6` anchor-isolation matrix certificate;
- independent C++ audits of the `(10,4)`, `(11,3)`, `(12,2)`, and `(13,1)`
  triangulation families;
- all 1,313,024 topology orbits and the exact partition
  `1,305,792 + 7,231 + 1`;
- the 913-node candidate tree and 16,075-node noncandidate forest;
- exact cross-stage root, graph, candidate-map, and residual-set links;
- rejection of corrupted roots, forces, orbit lists, malformed trees, and
  assertion-disabled execution.

No unresolved leaf or certificate inconsistency was found.

## Written-proof review

The Chinese proof explicitly supplies the geometric and topological links that
were important in earlier covering audits: ordinary Voronoi perturbation under
a strict radius buffer, irreducibility and private points, the boundary
no-return lemma, simplicity of the dual disk triangulation, stellar padding for
fewer than fourteen essential disks, the subgraph link used by local
isolation, and the topological covering argument for the two collapsed upper
faces.

This audit found no evident statement that must block repository publication.
That conclusion means the included exact replay succeeded and no concrete flaw
was identified in this review. It is not a claim of Lean/Coq formalization,
independent reimplementation, or completed external peer review.
