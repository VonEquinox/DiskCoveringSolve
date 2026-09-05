#!/usr/bin/env python3
"""Check the immutable Cover16 archive and regenerate its census in a temporary copy."""
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


def verify_manifest():
    manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    if len(manifest) != 90:
        raise RuntimeError("Expected 90 archived files")
    for rel, expected in manifest.items():
        name = Path(rel)
        if name.is_absolute() or ".." in name.parts or name.parts[0] != "proof_bundle":
            raise RuntimeError(f"Invalid manifest path: {rel}")
        digest = hashlib.sha256()
        with (ROOT / name).open("rb") as source:
            for block in iter(lambda: source.read(1 << 20), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    if (ROOT / "docs/r16_proof_zh.md").read_bytes() != (ROOT / "proof_bundle/PROOF_zh.md").read_bytes():
        raise RuntimeError("Standalone manuscript differs from archived proof")
    print("COVER16 ARCHIVE MANIFEST PASSED (90 files)", flush=True)


def environment():
    env = os.environ.copy()
    if sys.platform == "darwin" and shutil.which("brew"):
        result = subprocess.run(["brew", "--prefix", "boost"], capture_output=True, text=True)
        if result.returncode == 0:
            include = Path(result.stdout.strip()) / "include"
            if (include / "boost/multiprecision/cpp_int.hpp").is_file():
                env["CPLUS_INCLUDE_PATH"] = str(include) + (
                    os.pathsep + env["CPLUS_INCLUDE_PATH"] if env.get("CPLUS_INCLUDE_PATH") else ""
                )
    return env


def replay(jobs, report_dir=None):
    if report_dir is not None:
        report_dir = Path(report_dir).resolve()
        if report_dir.exists() or report_dir.is_relative_to((ROOT / "proof_bundle").resolve()):
            raise ValueError("Report directory must be new and outside the archived bundle")
    with tempfile.TemporaryDirectory(prefix="cover16_replay_") as tmp:
        work = Path(tmp) / "bundle"
        shutil.copytree(ROOT / "proof_bundle", work, ignore=shutil.ignore_patterns("__pycache__", "_runtime"))
        master = work / "reports/MASTER_VERIFIED.json"
        master.unlink(missing_ok=True)
        subprocess.run(
            [sys.executable, "-S", "-B", "verify_all.py", "--jobs", str(jobs)],
            cwd=work, env=environment(), check=True,
        )
        result = json.loads(master.read_text())
        if result.get("verified") is not True or result.get("unresolved_leaves") != 0:
            raise RuntimeError("Missing successful complete replay result")
        if result.get("enumeration_mode") != "regenerated_and_independently_audited":
            raise RuntimeError("Full census regeneration is required")
        verify_manifest()
        if report_dir is not None:
            shutil.copytree(work / "reports", report_dir)
    print("COVER16 EXACT MASTER REPLAY PASSED", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--jobs", type=int, choices=range(1, 6), default=1,
                        help="Concurrent census families; default 1 limits memory use")
    parser.add_argument("--report-dir", type=Path,
                        help="Optional new directory for fresh reports; runtime tables are discarded")
    args = parser.parse_args()
    verify_manifest()
    if not args.integrity_only:
        if sys.version_info < (3, 10):
            raise RuntimeError("Python 3.10 or newer is required; set PYTHON_BIN to that interpreter")
        replay(args.jobs, args.report_dir)


if __name__ == "__main__":
    main()
