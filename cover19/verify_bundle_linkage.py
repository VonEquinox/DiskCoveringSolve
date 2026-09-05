#!/usr/bin/env python3
"""Compare exact Cover19 geometry, stressed subgraphs and all survivor orbits."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import shutil

ROOT = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    folders = [ROOT / "proof_bundle", ROOT / "alternative_bundle"]
    sys.path.insert(0, str(folders[1] / "core"))
    import verify_peeling19 as peeling
    primary, alternative = [load(f"linked_geometry19_{i}", p / "core/geometry19.py")
                            for i, p in enumerate(folders)]
    assert primary.N == alternative.N == 19
    assert primary.B == alternative.B == 12
    assert primary.T == alternative.T == F(1, 13)
    assert primary.Q == alternative.ANCHORS
    assert len(set(primary.C)) == len(set(alternative.CENTERS)) == 19
    centers = [alternative.CENTERS.index(c) for c in primary.C]
    assert sorted(centers) == list(range(19)) and centers[:12] == list(range(12))
    full = {tuple(sorted(centers[v] for v in f)) for f in primary.FACES}
    active = {tuple(sorted(centers[v] for v in f)) for f in primary.ACTIVE}
    assert full == set(alternative.FACES) and active == set(alternative.ACTIVE)
    # The primary uses full-graph IDs; the alternative compresses stressed IDs.
    mapping = {i: i for i in range(12)}
    for i, j in enumerate(centers):
        if j in alternative.CENTER_MAP:
            mapping[12+i] = alternative.CENTER_MAP[j]
    for f in primary.ACTIVE:
        target = alternative.ACTIVE.index(tuple(sorted(centers[v] for v in f)))
        mapping[31+primary.FACES.index(f)] = 30+target
    assert sorted(mapping) == primary.ACTIVE_VERTICES
    assert sorted(mapping.values()) == list(range(36))
    assert all(primary.Z[i] == alternative.Z[j] for i, j in mapping.items())
    lookup = {frozenset(edge): w for edge, w in zip(alternative.EDGES, alternative.WEIGHTS)}
    mapped = {frozenset((mapping[u], mapping[v])): w for (u, v), w in zip(primary.EDGES, primary.WEIGHTS)}
    assert len(mapped) == len(lookup) == 42 and mapped == lookup
    assert primary.MU == F(-1, 156) and alternative.MUS == [primary.MU]*11
    omitted = [i for i in range(19) if 12+i not in mapping]
    assert omitted == [15] and centers[15] == 12 and primary.C[15] == (0, 0)
    # Explicitly bind the candidate to the tighter convex domain used only by
    # the primary lower bound: short consecutive gaps <= delta, cos(delta)=11/13.
    dots = [primary.dot(primary.Q[i], primary.Q[(i+1) % 12]) for i in range(12)]
    assert min(dots) >= F(11, 13)
    assert all(primary.det(primary.Q[i], primary.Q[(i+1) % 12]) > 0 for i in range(12))
    for folder, checker in zip(folders, ("verify_candidate19.py", "verify_upper19.py")):
        with tempfile.TemporaryDirectory(prefix="cover19_geometry_link_") as td:
            work = Path(td)
            shutil.copytree(folder / "core", work / "core", ignore=shutil.ignore_patterns("__pycache__", "*_verified.json"))
            subprocess.run([sys.executable, "-S", "-B", str(work / "core" / checker)],
                           cwd=work, check=True, capture_output=True, text=True)
    print("COVER19 EXACT GEOMETRY AND STRESS LINKAGE PASSED", flush=True)
    row_paths = [folders[0] / "enumeration/metric19_residuals.json",
                 folders[1] / "enumeration/survivors19.json"]
    rows = [json.loads(p.read_text()) for p in row_paths]
    families = {}
    for b in range(12, 19):
        codes = []
        for group in rows:
            selected = [r for r in group if r["B"] == b]
            canonical = dict(peeling.bfs_code(b, sorted(tuple(sorted(f)) for f in r["faces"]))
                             for r in selected)
            assert len(canonical) == len(selected), "Duplicate survivor orbit"
            codes.append(canonical)
        assert codes[0] == codes[1], (b, "Survivor sets or stabilizers differ")
        families[str(b)] = len(codes[0])
        print(f"COVER19 CROSS-BUNDLE ORBITS B={b}: {len(codes[0])} PASSED", flush=True)
    assert sum(families.values()) == 24127
    return {"verified": True, "exact_radius_squared": "1/13", "field": "Q(sqrt(3))",
            "geometry_checks_replayed": True, "center_map_primary_to_alternative": centers,
            "stressed_node_map": mapping, "stressed_rods": 42,
            "omitted_centers": {"primary": 15, "alternative": 12},
            "minimum_consecutive_anchor_dot": str(min(dots)),
            "candidate_in_primary_tight_convex_domain": True,
            "geometry_sha256": [digest(p / "core/geometry19.py") for p in folders],
            "survivor_files_sha256": [digest(p) for p in row_paths],
            "survivor_orbits_by_B": families,
            "scope": "Exact cross-bundle linkage supplements, not replaces, both full proofs."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is not None:
        target = args.output.resolve()
        if target.exists() or any(target.is_relative_to(ROOT / folder) for folder in ("proof_bundle", "alternative_bundle")):
            raise ValueError("Output must be new and outside immutable bundles")
    result = verify()
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("COVER19 EXACT GEOMETRY AND SURVIVOR LINKAGE PASSED", flush=True)
