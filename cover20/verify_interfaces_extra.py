"""Independently rebuild the Cover20 equations and bind the ZIP proof to them."""
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import re

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def face_blocks(text):
    blocks = re.findall(r"```text\n(.*?)\n```", text, re.S)
    return [ast.literal_eval("[" + block.strip().removesuffix(".") + "]") for block in blocks
            if block.strip().startswith("(0,")]


def verify(bundle):
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    bundle = Path(bundle)
    spec = importlib.util.spec_from_file_location("review_system20", bundle / "core/system20.py")
    system = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(system)
    data = json.loads((bundle / "core/root20.json").read_text())
    z = [F(int(v), int(data["Qx"])) for v in data["xnum"]]
    points = [[F(1), F(0)]] + [z[2*i:2*i+2] for i in range(51)]
    expected = [F(0)] * 198
    expected[95+102] = F(1)
    jac = [[F(0)] * 198 for _ in range(198)]
    # Reconstruct residuals and derivatives directly from squared lengths,
    # without using system20.CONS or its sparse-Hessian assembly.
    for row, (u, v) in enumerate(system.EDGES):
        diff = [points[u][d] - points[v][d] for d in range(2)]
        expected[row] = sum(x*x for x in diff) - z[102]
        jac[row][102] = -1
        weight = z[103+row]
        expected[95+102] -= weight
        jac[95+102][103+row] = -1
        for node, sign in ((u, 1), (v, -1)):
            if node == 0:
                continue
            for d in range(2):
                index = 2*(node-1)+d
                gradient = 2*sign*diff[d]
                jac[row][index] = gradient
                expected[95+index] += weight*gradient
                jac[95+index][103+row] = gradient
                for other, other_sign in ((u, 1), (v, -1)):
                    if other:
                        jac[95+index][2*(other-1)+d] += 2*sign*other_sign*weight
    for i in range(1, 13):
        row = 83+i-1
        expected[row] = sum(x*x for x in points[i]) - 1
        mu = z[103+row]
        for d in range(2):
            index = 2*(i-1)+d
            gradient = 2*points[i][d]
            jac[row][index] = gradient
            expected[95+index] += mu*gradient
            jac[95+index][103+row] = gradient
            jac[95+index][index] += 2*mu
    actual, actual_jac = system.evaluate(z, F(0))
    assert actual == expected
    assert actual_jac == jac

    inside = (bundle / "PROOF_zh.md").read_text()
    active, full = face_blocks(inside)
    assert active == system.ACTIVE and full == system.FACES
    original_hash = digest(bundle / "core/root20.json")
    inside_hashes = re.findall(r"(?m)^[0-9a-f]{64}$", inside)
    assert original_hash in inside_hashes
    rows = re.findall(r"(?m)^\| (\d+) \| ([+-]?[0-9.]+) \| ([+-]?[0-9.]+) \|$", inside)
    assert len(rows) == 20
    for label, x, y in rows:
        center = points[13+int(label)]
        assert max(abs(center[d]-F(value)) for d, value in enumerate((x, y))) < F(1, 10**15)
    report = {
        "verified": True,
        "independent_residual_components": 198,
        "independent_jacobian_entries": 198*198,
        "inside_proof_active_faces_match": True,
        "inside_proof_full_faces_match": True,
        "inside_display_centers_checked": 20,
        "actual_root_sha256": original_hash,
        "bundle_proof_sha256": digest(bundle / "PROOF_zh.md"),
        "scope": "Additional equation/document linkage, not a substitute for full replay or external review.",
    }
    print("COVER20 INDEPENDENT EQUATIONS AND ZIP PROOF LINKAGE PASSED", flush=True)
    return report


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent / "proof_bundle")
