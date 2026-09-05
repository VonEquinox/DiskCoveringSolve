# Cover20 release audit

Review and release preparation date: September 5, 2026. The maintainer chose
the ZIP's own proof as the published version. No blocking inconsistency was
found within that version in the completed direct review and full replay.
This is local repository acceptance, not independent external peer review.

## Source and version selection

The sole published proof source is `r20_verified_proof.zip`:

```text
ZIP bytes: 203343043
ZIP SHA-256:
2a3f92c971ddaa53da83188bfc1de5497bd3edda9e8e1cf7cd7fa5cb25301a1a
Root SHA-256:
9052d59645e7e888d64c8666c184d77d7e294c1e0b9a4b8914965b542c3e9ab9
```

All 51 original logical files are preserved. The original 37-entry accepting
manifest passed before and after direct replay. The ZIP contains no submitted
success reports, so all release reports were freshly produced. The repository
manifest also covers optional discovery sources and the original checksum file.

The separately submitted `PROOF_zh (1) 3.md` was reviewed but not selected.
It names a different root hash, uses `gamma=300`, `lambda=9/10000`, local radius
`13/100`, a 5,153-node candidate tree and a different integer force format.
The ZIP instead uses `gamma=200`, `lambda=11/25000`, radius `11/100` and a
6,171-node candidate tree. Its code and internal text agree. The unselected
document is not published or presented as a second accepted proof; matching
displayed radii alone would not establish exact correspondence of the versions.
Historical review metadata may record its hash, not certify its contents.

## Actual replay evidence

The direct master invocation was:

```bash
COVER20_WORKERS=4 /opt/homebrew/bin/python3 -S -B verify_all.py
```

It completed with exit code 0 in **448.59 seconds**, at
`2026-09-05T08:37:11.458529+00:00`. All 27 stages ran freshly from source.
The new repository runner also requires all these stages after reconstructing
the original large input. Fresh records of its release replay are kept in
`verification/repository_replay/`, separately from the direct audit in
`verification/independent_replay/`. "Independent replay" denotes local fresh
execution, not a review by an independent external party.

The repository reconstruction/full replay also completed successfully in
**537.29 seconds**, at `2026-09-05T08:56:12.547838+00:00`, with all 27 original
stages and the extra equation/text checker passing. Both runs report identical
accepting-input hashes, exact decimal intervals and global case totals. All
51 restored files were additionally byte-compared with the original ZIP.
Repository-wide input manifests for Cover11 through Cover20 and all seventeen
Cover20 wrapper regressions passed. This release did not repeat the full
mathematical replays of Cover11 through Cover19.

| Check | Direct full replay result |
|---|---:|
| Root system | 198 variables, two nested contraction boxes |
| Box half-widths | 10^-60 and 10^-80 |
| Upper complex | 33 vertices, 83 edges, 51 faces |
| Boundary caps | 13 |
| Exact-equality / strict-slack center faces | 19 / 6 |
| Geometric isolation matrix | 102 coordinates |
| Anchor displacement radius | 11/100 |
| Isolation margin | 1199/25000000 > 0 |
| Exhaustive families | 8, independently searched twice |
| Rooted survivors | 14,747,883 |
| Surviving orbit classes | 562,169 |
| Noncandidate cases | 562,168 |
| Noncandidate nodes / strict force leaves | 657,498 / 609,833 |
| Candidate nodes | 6,171 |
| Candidate splits / force leaves / local leaves | 3,085 / 2,984 / 102 |
| Unresolved leaves | 0 |
| Unpruned count regressions | 33 |
| Rejected malformed or invalid cases | 43 |
| Integer-lattice vs Fraction propagation tests | 100 |

The two enumerators visit 366,987,509 and 221,088,921 partial states,
respectively. Agreement is checked on full canonical survivor sets, not just
counts or digest strings. The retained sets are:

| B | I | Rooted survivors | Orbit classes |
|---:|---:|---:|---:|
| 12 | 8 | 3,710,028 | 154,931 |
| 13 | 7 | 5,827,029 | 224,321 |
| 14 | 6 | 3,870,440 | 138,451 |
| 15 | 5 | 1,189,410 | 39,713 |
| 16 | 4 | 150,976 | 4,753 |
| 17 | 3 | 0 | 0 |
| 18 | 2 | 0 | 0 |
| 19 | 1 | 0 | 0 |

These are outputs after necessary-condition partial-state pruning, not all
unpruned triangulations. The written peeling induction and pruning-safety
argument are essential; count agreement is not a completeness proof.

## Claim-to-evidence review

| Claim | Written argument | Reviewed executable evidence |
|---|---|---|
| Precisely defined isolated root and radius | ZIP Section 2 | Rational contraction on both boxes; rebuilt integer quadratic system |
| Entire unit disk is covered | Section 3 | Root-box upper checks, disk complex, all positive orientations, 13 arcs/caps |
| Every strictly smaller cover enters an allowed family | Sections 4-5 | Written Voronoi/private-point/no-return/dual/padding arguments; exact angular thresholds |
| Pruning preserves every hypothetical counterexample | Section 5 | Existing-edge short paths, safe long-arc branch selection, exact negative-cycle cuts |
| Peeling exhausts every admissible completion | Section 6 | Two fresh searches, full set comparison, extra unpruned recurrence checks |
| Noncandidate continuous domains are excluded | Sections 7 and 9 | Full streamed forest, exact force reconstruction, every strict separation inequality |
| Candidate local leaves apply to entire boxes | Sections 8-9 | Rational matrix perturbation/congruence, exact chord bounds, complete candidate tree |
| Candidate exemption matches the stored topology | Section 9.2 | Explicit full face-set permutation and 83-rod embedding into the 101-rod full graph |
| ZIP text specifies the actual equations | Sections 2-3 | Additional independent residual/Jacobian and text-to-data checker |

### Geometry and topology

The review checked the radius buffer before generic center perturbation,
ordinary restricted Voronoi containment, private points after removing
redundant disks, and the short-arc no-return proof. These establish the
simple disk dual; star subdivision pads fewer than twenty essential disks
without claiming the duplicated coordinates form a generic Voronoi diagram.
The radius threshold gives `12 <= B <= 19` and at least one interior cell.

For completeness of the dihedral quotient, reversing boundary labels may be
accompanied by reflecting all geometric points, relabeling anchors by their
boundary edges, and rotating to restore `q0=(1,0)`. This preserves lengths
and restores positive angular order. It is the geometric interpretation of
the reverse-root normalization in Section 6.3, not an added numerical claim.

The upper cover has no degenerate triangles. Its verified combinatorial
disk and positively oriented convex boundary give coverage of the anchor
polygon by the winding/extension argument in Section 3.4. The thirteen
circular caps then complete the unit disk. No visual inspection or sampling
is used to infer continuous coverage.

### Arithmetic and local/global linkage

The force representation keeps non-tree cycle forces, nonnegative support
multipliers and support vectors. Tree-edge forces are reconstructed exactly;
all 58 node balances are then recomputed independently. Outward trigonometric
intervals and integer square-root upper bounds certify strict inequalities
over the entire incoming angular box. No optimizer tolerance is accepted.

The local theorem acts on 102 geometric coordinates, not 198 KKT variables
or thirteen independent gaps. Its 95-by-102 Jacobian includes 83 rod and 12
circle constraints. Exact matrix perturbation bounds and a rational
congruence prove positivity with a quantified buffer. The theorem requires
only anchor closeness, not an unverified assumption that free centers are
already nearby. The candidate tree covers the whole necessary angular domain;
its local result alone is never used as a global conclusion.

The original candidate map identifies family B=13, index 217384. The full
25-face set is checked exactly, and the active 52-node graph is explicitly
embedded into the 58-node necessary graph. Discarding the six unused face
witnesses is a relaxation, valid for a lower bound.

The added checker reconstructs all 198 residuals and 39,204 Jacobian entries
directly from squared distances, without reusing the original sparse Hessian
assembly. It also checks the ZIP's full and active face tables, root hash,
and twenty displayed center coordinates. These checks passed in the direct
review and are included in the repository acceptance wrapper.

## Storage and software checks

The 197,440,381-byte compressed noncandidate forest exceeds the ordinary
GitHub single-file limit. Three lossless byte parts represent it; no force
record or original byte is removed or recompressed. The wrapper checks each
part and the original concatenated hash before invoking the unchanged master.
Direct execution inside the stored `proof_bundle/` is intentionally not the
repository entry point; reconstruction happens in a disposable copy.

Seventeen synthetic wrapper regressions cover reconstruction, corrupted and
missing parts, reordered and duplicate parts, wrong lengths, unsafe paths,
unexpected files, unsplit-input corruption, wrong proof versions, mutation
during restoration, incomplete stage sets, unresolved results, missing tests,
protected report paths, missing new success reports and export cleanup.
They supplement, not replace, the original 43 negative tests.

## Limits

Written convex geometry, planar topology, finite-dimensional analysis and
exhaustive-recursion induction remain trusted mathematics. Python exact
arithmetic and the C++ compiler/runtime remain software dependencies. The
two searches share the same underlying proof idea and are not external peer
reviews. No Lean/Coq formalization, historical priority, or classification of
all optimal covers is claimed. Optional numerical discovery tools are
preserved for methodology, but were not rerun for acceptance.
