PYTHON_BIN ?= python3

.PHONY: integrity quick full

integrity:
	cd proof_bundle && sha256sum -c SHA256SUMS

quick:
	PYTHON_BIN="$(PYTHON_BIN)" ./run_verification_portable.sh quick

full:
	PYTHON_BIN="$(PYTHON_BIN)" ./run_verification_portable.sh full
