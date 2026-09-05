# Archived, untrusted discovery sources

These files document how candidate and certificate data were proposed. They are **not imported by the proof verifier**. The accepting proof needs no numerical optimizer and no regeneration.

To experiment with regeneration, create a separate workspace. Copy these Python files, the candidate JSON, and the verifier Python modules from the parent directory into that workspace's top level. The scripts use their own directory as the working-data root. Do not run them against the only copy of the verified package: several generators intentionally overwrite data files.

The discovery environment is recorded in the parent `ENVIRONMENT.json`. In particular, `force_solver.py` uses SciPy SLSQP's returned multipliers. The standalone C++ enumerator takes boundary count, internal count and output path; its output is only proposed data. The independent `audit.cpp` and the root-edge count in `verify_partition.py` establish enumeration completeness.

The root JSON and multipliers from numerical discovery are not treated as exact answers. Only the independently replayed rational root enclosure, exact upper construction, integer-force trees and matrix congruence certificate are accepted.
