# Covering the Unit Disk with 11 Equal Disks

This directory contains the computer-assisted proof claim for covering the
unit disk with eleven congruent disks.

## Layout

- `proof_bundle/`: exact certificates, triangulation data, and verifier code;
- `verification/`: fresh full replay log and aggregate result;
- `docs/`: detailed mathematical method and verifier architecture;
- `run_verification.sh`: portable quick and full replay entry point.

## Replay

From the repository root:

```bash
make cover11-quick PYTHON_BIN="$(command -v python3)"
make cover11-full PYTHON_BIN="$(command -v python3)"
```

The full replay ends with:

```text
COVER11 FULL LEAF-BY-LEAF REPLAY PASSED
```
