#!/usr/bin/env python3
"""Recheck both Cover18 roots and compare their exact roots and survivor sets."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import permutations
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

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
    sys.path.insert(0, str(ROOT / "alternative_bundle/core"))
    import interval18 as iv
    import verify_peeling18 as peeling

    folders = [ROOT / "proof_bundle", ROOT / "alternative_bundle"]
    primary, alternative = [load(f"linked_system18_{i}", p / "core/system18.py")
                            for i, p in enumerate(folders)]
    for s in (primary, alternative):
        assert (s.B, s.N, s.G, s.P, s.M, s.D) == (12, 18, 90, 91, 72, 174)
    mappings = []
    for sign in (1, -1):
        for shift in range(12):
            for interior in permutations(range(12, 18)):
                centers = [(shift + sign*i) % 12 for i in range(12)] + list(interior)
                if {tuple(sorted(centers[v] for v in f)) for f in primary.FACES} != set(alternative.FACES):
                    continue
                if {tuple(sorted(centers[v] for v in f)) for f in primary.ACTIVE} == set(alternative.ACTIVE):
                    mappings.append((sign, shift, centers))
    assert mappings, "Full and active graph isomorphism missing"
    for folder in folders:
        with tempfile.TemporaryDirectory(prefix="cover18_root_link_") as tmp:
            for name in ("system18.py", "root18.json", "verify_root18.py"):
                shutil.copy2(folder / "core" / name, Path(tmp) / name)
            subprocess.run([sys.executable, "-S", "-B", "verify_root18.py"], cwd=tmp,
                           check=True, capture_output=True, text=True)
    paths = [p / "core/root18.json" for p in folders]
    raw = [json.loads(p.read_text()) for p in paths]
    middle, other = [[F(int(v), int(d["Qx"])) for v in d["xnum"]] for d in raw]
    eps = F(1, 10**80)
    box = [(v-eps, v+eps) for v in other]
    points = [[iv.point(1), iv.point(0)]] + [box[2*i:2*i+2] for i in range(45)]
    # Restore the omitted q0 multiplier; total torque makes its force radial.
    mu0 = iv.point(0)
    for k, (u, v) in enumerate(alternative.EDGES):
        if 0 in (u, v):
            dx = iv.sub(iv.point(1), points[v if u == 0 else u][0])
            mu0 = iv.sub(mu0, iv.mul(box[91+k], dx))
    multipliers = [mu0] + box[163:]
    transports = []
    for sign, shift, centers in mappings:
        anchors = [(shift + sign*i - (sign == -1)) % 12 for i in range(12)]
        witnesses = [alternative.ACTIVE.index(tuple(sorted(centers[v] for v in f)))
                     for f in primary.ACTIVE]
        nodes = anchors + [12+i for i in centers] + [30+i for i in witnesses]
        assert sorted(nodes) == list(range(46))
        lookup = {frozenset(e): k for k, e in enumerate(alternative.EDGES)}
        rods = [lookup[frozenset((nodes[u], nodes[v]))] for u, v in primary.EDGES]
        assert sorted(rods) == list(range(72))
        anchor = points[anchors[0]]
        image = []
        for node in nodes[1:]:
            image += [iv.dot(anchor, points[node]), iv.scale(iv.cross(anchor, points[node]), sign)]
        image += [box[90]] + [box[91+k] for k in rods]
        image += [multipliers[anchors[i]] for i in range(1, 12)]
        assert len(image) == 174
        error = max(max(abs(lo-c), abs(hi-c)) for (lo, hi), c in zip(image, middle))
        if error < F(1, 10**60):
            transports.append({"sign": sign, "shift": shift, "center_map": centers,
                               "anchor_map": anchors, "witness_map": witnesses,
                               "maximum_transport_deviation": str(error),
                               "uniqueness_radius": str(F(1, 10**60))})
    assert transports, "Exact root transport misses the primary uniqueness box"
    print("COVER18 EXACT ROOT TRANSPORT PASSED", flush=True)
    row_paths = [folders[0] / "enumeration/metric18_residuals.json",
                 folders[1] / "enumeration/survivors18.json"]
    rows = [json.loads(p.read_text()) for p in row_paths]
    families = {}
    for b in range(11, 18):
        codes = []
        for group in rows:
            selected = [r for r in group if r["B"] == b]
            canonical = dict(peeling.bfs_code(b, sorted(tuple(sorted(f)) for f in r["faces"]))
                             for r in selected)
            assert len(canonical) == len(selected), "Duplicate survivor orbit"
            codes.append(canonical)
        assert codes[0] == codes[1], (b, "Survivor orbit sets or stabilizers differ")
        families[str(b)] = len(codes[0])
        print(f"COVER18 CROSS-BUNDLE ORBITS B={b}: {len(codes[0])} PASSED", flush=True)
    assert sum(families.values()) == 16220
    return {"verified": True, "both_root_certificates_replayed": True,
            "same_exact_t_after_orthogonal_transport": True,
            "root_sha256": [digest(p) for p in paths], "valid_root_transports": transports,
            "survivor_orbits_by_B": families, "survivor_files_sha256": [digest(p) for p in row_paths],
            "scope": "Additional exact linkage, not a replacement for either full replay."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is not None:
        target = args.output.resolve()
        if target.exists() or any(target.is_relative_to(ROOT / folder) for folder in ("proof_bundle", "alternative_bundle")):
            raise ValueError("Output must be new and outside the immutable bundles")
    result = verify()
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("COVER18 EXACT ROOT AND SURVIVOR LINKAGE PASSED", flush=True)
