#!/usr/bin/env python3
"""Extra unpruned peeling regressions with four to six interior vertices."""
import argparse
import json
from math import factorial
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "proof_bundle/core"))
from verify_enumeration17 import labelled_count


def run(command):
    return subprocess.run(list(map(str, command)), check=True, capture_output=True, text=True).stdout


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sanitize", action="store_true", help="Enable AddressSanitizer and UndefinedBehaviorSanitizer")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("A C++17 compiler is required")
    start = time.time()
    results = []
    with tempfile.TemporaryDirectory(prefix="cover17_extra_counts_") as td:
        tmp = Path(td)
        options = ["-std=c++17", "-O1"]
        if args.sanitize:
            options += ["-fsanitize=address,undefined", "-fno-sanitize-recover=all", "-fno-omit-frame-pointer"]
        for tool in ["peel_enum", "peel_audit"]:
            run([compiler, *options, ROOT / "proof_bundle/enumeration" / (tool + ".cpp"), "-o", tmp / tool])
        for b in [3, 4]:
            for i in [4, 5, 6]:
                n = b + i
                prefix = tmp / f"N{n}_B{b}"
                report = tmp / f"audit_N{n}_B{b}.json"
                run([tmp / "peel_enum", n, b, 9650, 20367, 35325, prefix, 0])
                run([tmp / "peel_audit", n, b, 9650, 20367, 35325, prefix.with_suffix(".txt"), report, 0])
                primary = json.loads(prefix.with_suffix(".json").read_text())
                audit = json.loads(report.read_text())
                count = labelled_count(b, i)
                assert primary["complete"] is True and primary["pruning"] is False
                assert audit["verified"] is True
                assert primary["rooted_survivors"] * factorial(i) == count
                assert audit["rooted_survivors"] == primary["rooted_survivors"]
                assert audit["surviving_orbits"] == primary["orbits"]
                results.append({"N": n, "B": b, "I": i, "labelled_count": count,
                                "rooted_survivors": primary["rooted_survivors"], "orbits": primary["orbits"]})
                print(f"EXTRA UNPRUNED COUNT PASSED N={n} B={b} I={i}", flush=True)
    result = {"verified": True, "sanitizers": args.sanitize, "cases": results, "seconds": time.time() - start}
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + "\n")
    print("COVER17 EXTRA UNPRUNED COUNTS PASSED", flush=True)


if __name__ == "__main__":
    main()
