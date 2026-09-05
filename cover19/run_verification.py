#!/usr/bin/env python3
"""Check immutable Cover19 inputs and replay both chains in disposable copies."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
BUNDLES = {"primary": "proof_bundle", "alternative": "alternative_bundle"}
COUNTS = {"primary": 23, "alternative": 66}
STEPS = {
    "primary": ["verify_interfaces19.py", "verify_candidate19.py",
                "verify_enumeration19.py", "verify_forest19.py",
                "test_fail_closed19.py"],
    "alternative": ["verify_interfaces19", "verify_upper19", "anchor_isolation19",
                    "compile_dfs", "compile_bfs"]
                   + [f"{mode}_B{b}" for b in range(12, 19) for mode in ("dfs", "bfs")]
                   + ["verify_peeling19", "verify_forest19", "verify_recursion_tests19", "test_fail_closed19"],
}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(stream.read()).hexdigest()


def verify_manifest():
    manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    if len(manifest) != sum(COUNTS.values()):
        raise RuntimeError("Expected 89 archived files")
    for rel, expected in manifest.items():
        name = Path(rel)
        if name.is_absolute() or ".." in name.parts or not name.parts or name.parts[0] not in BUNDLES.values():
            raise RuntimeError(f"Invalid manifest path: {rel}")
        if digest(ROOT / name) != expected:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    for variant, folder in BUNDLES.items():
        actual = {str(p.relative_to(ROOT)) for p in (ROOT / folder).rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts}
        expected = {p for p in manifest if Path(p).parts[0] == folder}
        if actual != expected or len(expected) != COUNTS[variant]:
            raise RuntimeError(f"Unexpected archive file set: {variant}")
        proof = "r19_proof_zh.md" if variant == "primary" else "r19_alternative_proof_zh.md"
        if (ROOT / "docs" / proof).read_bytes() != (ROOT / folder / "PROOF_zh.md").read_bytes():
            raise RuntimeError(f"Standalone proof differs: {variant}")
    print("COVER19 ARCHIVE MANIFEST PASSED (89 files)", flush=True)


def validate_result(variant, result):
    if result.get("verified") is not True or result.get("unresolved_leaves") != 0:
        raise RuntimeError("Missing successful complete replay result")
    key = "step" if variant == "primary" else "stage"
    if [s.get(key) for s in result.get("steps", [])] != STEPS[variant]:
        raise RuntimeError("Every accepting stage must run in order")
    orbit_key = "surviving_topology_orbits" if variant == "primary" else "surviving_orbits"
    radius_key = "exact_squared_radius" if variant == "primary" else "exact_radius_squared"
    if result.get(radius_key) != "1/13" or result.get(orbit_key) != 24127:
        raise RuntimeError("Unexpected exact radius or surviving-set size")
    case_key = "noncandidate_cases" if variant == "primary" else "cases"
    if result.get("noncandidate_forest", {}).get(case_key) != 24126:
        raise RuntimeError("Missing noncandidate cases")
    if variant == "primary":
        negatives = result.get("rejection_tests", {}).get("rejection_tests")
        unpruned = result.get("unpruned_recurrence_tests")
        if result.get("candidate_global_lower", {}).get("positive_rods") != 42 or result.get("candidate_cases") != 1:
            raise RuntimeError("Missing global candidate lower bound")
    else:
        negatives = result.get("negative_tests", {}).get("rejection_tests")
        unpruned = result.get("unpruned_recurrence_tests", {}).get("cases")
        if result.get("positive_stress_graph", {}).get("geometric_dimension") != 70:
            raise RuntimeError("Missing stressed-graph isolation")
    if negatives != (35 if variant == "primary" else 36) or unpruned != (26 if variant == "primary" else 29):
        raise RuntimeError("Missing regression tests")


def report_target(report_dir):
    if report_dir is None:
        return None
    target = Path(report_dir).resolve()
    if target.exists() or any(target.is_relative_to((ROOT / folder).resolve()) for folder in BUNDLES.values()):
        raise ValueError("Report directory must be new and outside both archived bundles")
    return target


def replay(variant, report_dir=None):
    report_dir = report_target(report_dir)
    with tempfile.TemporaryDirectory(prefix=f"cover19_{variant}_") as td:
        work = Path(td) / "bundle"
        shutil.copytree(ROOT / BUNDLES[variant], work, ignore=shutil.ignore_patterns("__pycache__", "_runtime"))
        if (work / "reports").exists():
            shutil.rmtree(work / "reports")
        (work / "reports").mkdir()
        for old in (work / "core").glob("*_verified.json"):
            old.unlink()
        env = dict(os.environ)
        env.pop("PYTHONOPTIMIZE", None)
        env.pop("PYTHONPATH", None)
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        subprocess.run([sys.executable, "-S", "-B", "verify_all.py"], cwd=work, check=True, env=env)
        result = json.loads((work / "reports/MASTER_VERIFIED.json").read_text())
        validate_result(variant, result)
        verify_manifest()
        if report_dir is not None:
            shutil.copytree(work / "reports", report_dir)
            (report_dir / "core").mkdir()
            for fresh in (work / "core").glob("*_verified.json"):
                shutil.copyfile(fresh, report_dir / "core" / fresh.name)
    print(f"COVER19 {variant.upper()} EXACT MASTER REPLAY PASSED", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=[*BUNDLES, "both"], default="both")
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--report-dir", type=Path, help="New directory for fresh reports and exact linkage")
    args = parser.parse_args()
    verify_manifest()
    if args.integrity_only:
        return
    if sys.version_info < (3, 10):
        raise RuntimeError("Python 3.10+ required")
    target = report_target(args.report_dir)
    variants = list(BUNDLES) if args.variant == "both" else [args.variant]
    for variant in variants:
        replay(variant, target / variant if target is not None else None)
    from verify_bundle_linkage import verify
    linkage = verify()
    verify_manifest()
    if target is not None:
        (target / "linkage_verified.json").write_text(json.dumps(linkage, indent=2) + "\n")
    print("COVER19 REQUESTED REPLAYS AND EXACT LINKAGE PASSED", flush=True)


if __name__ == "__main__":
    main()
