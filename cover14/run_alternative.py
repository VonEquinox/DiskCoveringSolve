#!/usr/bin/env python3
"""Replay the archived alternative package in a disposable working copy."""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT / "alternative_bundle"


def verify_manifest():
    seen = set()
    for line in (BUNDLE / "SHA256SUMS.txt").read_text().splitlines():
        expected, rel = line.split("  ", 1)
        name = Path(rel)
        if name.is_absolute() or ".." in name.parts or rel in seen:
            raise RuntimeError(f"Invalid or duplicate manifest path: {rel}")
        seen.add(rel)
        if hashlib.sha256((BUNDLE / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    if len(seen) != 64:
        raise RuntimeError("Expected 64 archived manifest entries")
    print("COVER14 ALTERNATIVE MANIFEST PASSED (64 files)", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    args = parser.parse_args()
    verify_manifest()
    if args.integrity_only:
        return
    env = os.environ.copy()
    if sys.platform == "darwin" and shutil.which("brew"):
        result = subprocess.run(["brew", "--prefix", "boost"], capture_output=True, text=True)
        if result.returncode == 0:
            include = Path(result.stdout.strip()) / "include"
            if (include / "boost/multiprecision/cpp_int.hpp").is_file():
                env["CPLUS_INCLUDE_PATH"] = str(include) + (
                    os.pathsep + env["CPLUS_INCLUDE_PATH"] if env.get("CPLUS_INCLUDE_PATH") else ""
                )
    with tempfile.TemporaryDirectory(prefix="cover14_alternative_") as tmp:
        work = Path(tmp) / "bundle"
        shutil.copytree(BUNDLE, work)
        subprocess.run([sys.executable, "-S", "-B", "verify_all.py"], cwd=work, env=env, check=True)
    subprocess.run([sys.executable, "-S", "-B", str(ROOT / "verify_bundle_linkage.py")], check=True)
    verify_manifest()
    print("COVER14 ALTERNATIVE MASTER REPLAY AND ROOT LINKAGE PASSED", flush=True)


if __name__ == "__main__":
    main()
