PYTHON_BIN ?= python3

.PHONY: integrity quick full cover11-integrity cover11-quick cover11-full cover12-integrity cover12-full cover13-integrity cover13-full

integrity: cover11-integrity cover12-integrity cover13-integrity

quick: cover11-quick

full: cover11-full cover12-full cover13-full

cover11-integrity:
	cd cover11/proof_bundle && sha256sum -c SHA256SUMS

cover11-quick:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover11/run_verification.sh quick

cover11-full:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover11/run_verification.sh full

cover12-integrity:
	cd cover12/proof_bundle && { command -v sha256sum >/dev/null 2>&1 && sha256sum -c SHA256SUMS || shasum -a 256 -c SHA256SUMS; }

cover12-full:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover12/run_verification.sh

cover13-integrity:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover13/verify_manifest.py

cover13-full:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover13/run_verification.sh
