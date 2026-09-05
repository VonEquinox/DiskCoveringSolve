PYTHON_BIN ?= python3

.PHONY: integrity quick full cover11-integrity cover11-quick cover11-full cover12-integrity cover12-full cover13-integrity cover13-full cover14-integrity cover14-full
.PHONY: cover14-alternative-integrity cover14-alternative-full
.PHONY: cover15-integrity cover15-full cover15-primary-full cover15-alternative-full
.PHONY: cover16-integrity cover16-full cover16-runner-tests
.PHONY: cover17-integrity cover17-full cover17-runner-tests cover17-extra-tests
.PHONY: cover18-integrity cover18-full cover18-primary-full cover18-alternative-full cover18-runner-tests cover18-linkage
.PHONY: cover19-integrity cover19-full cover19-primary-full cover19-alternative-full cover19-runner-tests cover19-linkage
.PHONY: cover20-integrity cover20-full cover20-runner-tests

integrity: cover11-integrity cover12-integrity cover13-integrity cover14-integrity cover14-alternative-integrity cover15-integrity cover16-integrity cover17-integrity cover18-integrity cover19-integrity cover20-integrity

quick: cover11-quick

full: cover11-full cover12-full cover13-full cover14-full cover14-alternative-full cover15-full cover16-full cover17-full cover18-full cover19-full cover20-full

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

cover17-integrity:
	"$(PYTHON_BIN)" -S -B cover17/run_verification.py --integrity-only

cover17-full:
	"$(PYTHON_BIN)" -S -B cover17/run_verification.py

cover17-runner-tests:
	"$(PYTHON_BIN)" -S -B cover17/test_run_verification.py

cover17-extra-tests:
	"$(PYTHON_BIN)" -S -B cover17/test_peeling_counts.py

cover18-integrity:
	"$(PYTHON_BIN)" -S -B cover18/run_verification.py --integrity-only

cover18-full:
	"$(PYTHON_BIN)" -S -B cover18/run_verification.py --variant both

cover18-primary-full:
	"$(PYTHON_BIN)" -S -B cover18/run_verification.py --variant primary

cover18-alternative-full:
	"$(PYTHON_BIN)" -S -B cover18/run_verification.py --variant alternative

cover18-runner-tests:
	"$(PYTHON_BIN)" -S -B cover18/test_run_verification.py

cover18-linkage:
	"$(PYTHON_BIN)" -S -B cover18/verify_bundle_linkage.py

cover19-integrity:
	"$(PYTHON_BIN)" -S -B cover19/run_verification.py --integrity-only

cover19-full:
	"$(PYTHON_BIN)" -S -B cover19/run_verification.py --variant both

cover19-primary-full:
	"$(PYTHON_BIN)" -S -B cover19/run_verification.py --variant primary

cover19-alternative-full:
	"$(PYTHON_BIN)" -S -B cover19/run_verification.py --variant alternative

cover19-runner-tests:
	"$(PYTHON_BIN)" -S -B cover19/test_run_verification.py

cover19-linkage:
	"$(PYTHON_BIN)" -S -B cover19/verify_bundle_linkage.py

cover20-integrity:
	"$(PYTHON_BIN)" -S -B cover20/run_verification.py --integrity-only

cover20-full:
	"$(PYTHON_BIN)" -S -B cover20/run_verification.py

cover20-runner-tests:
	"$(PYTHON_BIN)" -S -B cover20/test_run_verification.py
