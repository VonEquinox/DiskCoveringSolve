#!/usr/bin/env python3
"""Verify Cover17's immutable archive and replay every stage in a disposable copy."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
STEPS = ["verify_interfaces17.py", "verify_root17.py", "verify_upper17.py",
         "anchor_isolation17.py", "verify_enumeration17.py", "verify_forest17.py",
         "test_fail_closed17.py"]


def verify_manifest():
    manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    if len(manifest) != 46:
        raise RuntimeError("Expected 46 archived files")
    for rel, expected in manifest.items():
        name = Path(rel)
        if name.is_absolute() or ".." in name.parts or name.parts[0] != "proof_bundle":
            raise RuntimeError(f"Invalid manifest path: {rel}")
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    if (ROOT / "docs/r17_proof_zh.md").read_bytes() != (ROOT / "proof_bundle/PROOF_zh.md").read_bytes():
        raise RuntimeError("Standalone proof differs from archived proof")
    print("COVER17 ARCHIVE MANIFEST PASSED (46 files)", flush=True)


def replay(report_dir=None):
    if report_dir is not None:
        report_dir = Path(report_dir).resolve()
        if report_dir.exists() or report_dir.is_relative_to((ROOT / "proof_bundle").resolve()):
            raise ValueError("Report directory must be new and outside the archived bundle")
    with tempfile.TemporaryDirectory(prefix="cover17_replay_") as td:
        work = Path(td) / "bundle"
        shutil.copytree(ROOT / "proof_bundle", work, ignore=shutil.ignore_patterns("__pycache__"))
        master = work / "reports/MASTER_VERIFIED.json"
        master.unlink(missing_ok=True)
        for old in (work / "core").glob("*_verified.json"):
            old.unlink()
        subprocess.run([sys.executable, "-S", "-B", "verify_all.py"], cwd=work, check=True)
        result = json.loads(master.read_text())
        if result.get("verified") is not True or result.get("unresolved_leaves") != 0:
            raise RuntimeError("Missing successful complete replay result")
        if [step["step"] for step in result.get("steps", [])] != STEPS:
            raise RuntimeError("All seven accepting stages must run")
        if result.get("kkt_dimension") != 157 or result.get("surviving_topology_orbits") != 1344:
            raise RuntimeError("Unexpected system or surviving-set size")
        verify_manifest()
        if report_dir is not None:
            shutil.copytree(work / "reports", report_dir)
            (report_dir / "core").mkdir()
            for fresh in (work / "core").glob("*_verified.json"):
                shutil.copyfile(fresh, report_dir / "core" / fresh.name)
    print("COVER17 EXACT MASTER REPLAY PASSED", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--report-dir", type=Path, help="Optional new directory for fresh stage reports")
    args = parser.parse_args()
    verify_manifest()
    if not args.integrity_only:
        if sys.version_info < (3, 10):
            raise RuntimeError("Python 3.10+ required; set PYTHON_BIN to that interpreter")
        replay(args.report_dir)


if __name__ == "__main__":
    main()
