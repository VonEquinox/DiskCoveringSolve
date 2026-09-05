PYTHON_BIN ?= python3

.PHONY: integrity quick full cover11-integrity cover11-quick cover11-full cover12-integrity cover12-full cover13-integrity cover13-full cover14-integrity cover14-full
.PHONY: cover14-alternative-integrity cover14-alternative-full
.PHONY: cover15-integrity cover15-full cover15-primary-full cover15-alternative-full
.PHONY: cover16-integrity cover16-full cover16-runner-tests

integrity: cover11-integrity cover12-integrity cover13-integrity cover14-integrity cover14-alternative-integrity cover15-integrity cover16-integrity

quick: cover11-quick

full: cover11-full cover12-full cover13-full cover14-full cover14-alternative-full cover15-full cover16-full

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

cover14-integrity:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover14/verify_manifest.py

cover14-full:
	PYTHON_BIN="$(PYTHON_BIN)" ./cover14/run_verification.sh

cover14-alternative-integrity:
	"$(PYTHON_BIN)" -S -B cover14/run_alternative.py --integrity-only

cover14-alternative-full:
	"$(PYTHON_BIN)" -S -B cover14/run_alternative.py

cover15-integrity:
	"$(PYTHON_BIN)" -S -B cover15/run_verification.py --integrity-only

cover15-full:
	"$(PYTHON_BIN)" -S -B cover15/run_verification.py --variant both

cover15-primary-full:
	"$(PYTHON_BIN)" -S -B cover15/run_verification.py --variant primary

cover15-alternative-full:
	"$(PYTHON_BIN)" -S -B cover15/run_verification.py --variant alternative

cover16-integrity:
	"$(PYTHON_BIN)" -S -B cover16/run_verification.py --integrity-only

cover16-full:
	"$(PYTHON_BIN)" -S -B cover16/run_verification.py

cover16-runner-tests:
	"$(PYTHON_BIN)" -S -B cover16/test_run_verification.py
