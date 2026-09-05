#!/usr/bin/env python3
"""Recheck both roots and identify them by exact rotation and relabeling."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

if not __debug__:
    raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "proof_bundle/core"))
import system15 as primary
import interval15 as interval

spec = importlib.util.spec_from_file_location("alternative_system15", ROOT / "alternative_bundle/core/system15.py")
alternative = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alternative)


def upper_connected(system):
    b, n = system.B, system.N
    faces = [tuple(b + v for v in face) for face in system.FACES]
    for i in range(b):
        faces.extend([((i-1) % b, i, b+i), (i, b+(i+1) % b, b+i)])
    adjacency = {v: set() for v in range(b+n)}
    for face in faces:
        for u in face:
            adjacency[u].update(set(face)-{u})
    reached, todo = {0}, [0]
    while todo:
        for v in adjacency[todo.pop()]:
            if v not in reached:
                reached.add(v)
                todo.append(v)
    assert reached == set(adjacency), "Upper complex disconnected"
    return len(reached)


def transport():
    assert (primary.B, primary.N, primary.G, primary.P, primary.M, primary.D) == (11, 15, 76, 77, 61, 148)
    assert (alternative.B, alternative.N, alternative.G, alternative.P, alternative.M, alternative.D) == (11, 15, 76, 77, 61, 148)
    centers = [5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 13, 14, 11, 12]
    anchors = [(i+5) % 11 for i in range(11)]
    assert sorted(centers) == list(range(15))
    assert {tuple(sorted(centers[i] for i in face)) for face in primary.FACES} == set(alternative.FACES)
    active = {tuple(face): k for k, face in enumerate(alternative.ACTIVE)}
    witnesses = [active[tuple(sorted(centers[i] for i in face))] for face in primary.ACTIVE]
    nodes = anchors + [11+i for i in centers] + [26+i for i in witnesses]
    assert sorted(nodes) == list(range(39))
    lookup = {frozenset(edge): k for k, edge in enumerate(alternative.EDGES)}
    rods = [lookup[frozenset((nodes[u], nodes[v]))] for u, v in primary.EDGES]
    assert sorted(rods) == list(range(61))

    paths = [ROOT / folder / "core/root15.json" for folder in ["proof_bundle", "alternative_bundle"]]
    old, new = [json.loads(p.read_text()) for p in paths]
    middle = [F(int(v), int(old["Qx"])) for v in old["xnum"]]
    other = [F(int(v), int(new["Qx"])) for v in new["xnum"]]
    assert len(middle) == len(other) == 148
    radius = F(1, 10**80)
    box = [(v-radius, v+radius) for v in other]
    points = [[interval.point(1), interval.point(0)]] + [box[2*i:2*i+2] for i in range(38)]
    # Recover the omitted anchor multiplier. All other node forces are radial
    # or zero; cancellation of total torque makes this remaining force radial.
    mu0 = interval.point(0)
    for k, (u, v) in enumerate(alternative.EDGES):
        if 0 in (u, v):
            displacement = interval.sub(interval.point(1), points[v if u == 0 else u][0])
            mu0 = interval.sub(mu0, interval.mul(box[77+k], displacement))
    multipliers = [mu0] + box[138:]
    anchor = points[anchors[0]]
    image = []
    for node in nodes[1:]:
        image += [interval.dot(anchor, points[node]), interval.cross(anchor, points[node])]
    image += [box[76]] + [box[77+k] for k in rods]
    image += [multipliers[anchors[i]] for i in range(1, 11)]
    assert len(image) == 148
    error = max(max(abs(lo-c), abs(hi-c)) for (lo, hi), c in zip(image, middle))
    assert error < F(1, 10**60), "Root transport misses primary uniqueness box"
    return {
        "verified": True,
        "center_map_primary_to_alternative": centers,
        "anchor_map_primary_to_alternative": anchors,
        "witness_map_primary_to_alternative": witnesses,
        "maximum_transport_deviation": str(error),
        "primary_uniqueness_radius": str(F(1, 10**60)),
        "primary_root_sha256": hashlib.sha256(paths[0].read_bytes()).hexdigest(),
        "alternative_root_sha256": hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        "connected_upper_vertices": [upper_connected(s) for s in [primary, alternative]],
    }


def verify():
    for folder in ["proof_bundle", "alternative_bundle"]:
        with tempfile.TemporaryDirectory(prefix="cover15_root_") as tmp:
            work = Path(tmp)
            for name in ["root15.json", "system15.py", "verify_root15.py"]:
                shutil.copy2(ROOT / folder / "core" / name, work / name)
            result = subprocess.run([sys.executable, "-S", "-B", "verify_root15.py"], cwd=work, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(f"{folder} root replay failed:\n{result.stdout}\n{result.stderr}")
    report = transport()
    report["both_root_certificates_replayed"] = True
    report["conclusion"] = "Both isolated roots have the same t coordinate after graph relabeling and rotation."
    return report


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
