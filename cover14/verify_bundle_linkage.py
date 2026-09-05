#!/usr/bin/env python3
"""Transport the alternative KKT root into the primary uniqueness box.

Root existence and uniqueness are established by the two full replays. This
check binds their graph labels and rational boxes, including the multiplier of
the anchor which becomes fixed after rotation.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "proof_bundle/core"))
sys.path.insert(0, str(ROOT / "alternative_bundle"))
import system14 as primary
import framework as alternative
from intervals import I, dot, pt


def verify():
    centers = [3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 11, 12, 13, 10]
    anchors = [(i + 3) % 10 for i in range(10)]
    assert sorted(centers) == list(range(14))
    assert {tuple(sorted(centers[i] for i in f)) for f in primary.FACES} == set(alternative.ALL_FACES)
    face_lookup = {tuple(f): k for k, f in enumerate(alternative.ACTIVE_FACES)}
    witnesses = [face_lookup[tuple(sorted(centers[i] for i in f))] for f in primary.ACTIVE]
    nodes = anchors + [10 + i for i in centers] + [24 + i for i in witnesses]
    assert sorted(nodes) == list(range(32))
    edge_lookup = {frozenset(e): k for k, e in enumerate(alternative.EDGES)}
    edges = [edge_lookup[frozenset((nodes[u], nodes[v]))] for u, v in primary.EDGES]
    assert sorted(edges) == list(range(44))
    assert primary.TID == alternative.TID == 62

    old_path = ROOT / "proof_bundle/core/root14.json"
    new_path = ROOT / "alternative_bundle/root_certificate.json"
    old = json.loads(old_path.read_text())
    new = json.loads(new_path.read_text())
    middle = [F(int(v), int(old["Qx"])) for v in old["xnum"]]
    other = [F(int(v), int(new["Qx"])) for v in new["xnum"]]
    radius = F(int(new["rho_num"]), int(new["rho_den"]))
    assert radius > 0 and len(middle) == len(other) == 116
    box = [I(v - radius, v + radius) for v in other]
    points = [pt(1, 0)] + [(box[2*i], box[2*i+1]) for i in range(31)]

    # At a KKT root the omitted anchor force is radial: pairwise cancellation
    # gives zero total torque, and every other node is already stationary.
    mu0 = -sum(
        (box[63+k] * (1-points[v if u == 0 else u][0])
         for k, (u, v) in enumerate(alternative.EDGES) if 0 in (u, v)),
        I.point(0),
    )
    multipliers = [mu0] + box[107:]
    q = points[anchors[0]]

    def rotate(p):
        return dot(q, p), q[0]*p[1] - q[1]*p[0]

    # The selected anchor is exactly unit at the true root, so this is a
    # rotation taking that anchor to (1,0); all other points use the same map.
    image = [v for node in nodes[1:] for v in rotate(points[node])]
    image += [box[62]] + [box[63+k] for k in edges]
    image += [multipliers[anchors[i]] for i in range(1, 10)]
    assert len(image) == 116
    deviation = max(max(abs(iv.lo-c), abs(iv.hi-c)) for iv, c in zip(image, middle))
    uniqueness_radius = F(1, 10**60)
    assert deviation < uniqueness_radius, "Transport misses primary uniqueness box"
    return {
        "verified": True,
        "meaning": "After both root replays, rotation and relabeling identify the same KKT root and t coordinate.",
        "center_map_primary_to_alternative": centers,
        "anchor_map_primary_to_alternative": anchors,
        "witness_map_primary_to_alternative": witnesses,
        "maximum_transport_deviation": str(deviation),
        "primary_uniqueness_radius": str(uniqueness_radius),
        "primary_root_sha256": hashlib.sha256(old_path.read_bytes()).hexdigest(),
        "alternative_root_sha256": hashlib.sha256(new_path.read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
