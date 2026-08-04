PYTHON ?= python3

.PHONY: all extract persistence cluster robustness sensitivity thermal regression test clean

all: extract persistence cluster robustness sensitivity thermal regression

extract:
	$(PYTHON) scripts/01_extract_spectra.py

persistence:
	$(PYTHON) scripts/02_compute_persistence.py

cluster:
	$(PYTHON) scripts/03_clustering.py

robustness:
	$(PYTHON) scripts/04_robustness.py

sensitivity:
	$(PYTHON) scripts/05_sensitivity.py

thermal:
	$(PYTHON) scripts/06_thermal_targets.py

regression:
	$(PYTHON) scripts/07_property_regression.py

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf data/processed/* results/* figures/*
