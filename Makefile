PYTHON ?= python3

.PHONY: all extract persistence cluster robustness sensitivity thermal regression multimodal peaks alignment test clean

all: extract persistence cluster robustness sensitivity thermal regression multimodal peaks alignment

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

multimodal:
	$(PYTHON) scripts/08_multimodal.py

# The controlled comparison against conventional peak descriptors. Parts B and
# D are the expensive ones; they can be run a target or an artefact at a time
# with --targets / --artefacts, and each invocation adds its rows to the table.
peaks:
	$(PYTHON) scripts/09_peak_features.py
	$(PYTHON) scripts/09b_peak_figure.py

# Pre-processing and alignment baselines. The noise artefact is run twice, with
# the fixed and the adaptive pruning threshold, since Sect. 5 prescribes the
# adaptive one for exactly that case.
alignment:
	$(PYTHON) scripts/10_alignment.py --parts B --targets T5
	$(PYTHON) scripts/10_alignment.py --parts D --artefacts "wavenumber shift"
	$(PYTHON) scripts/10_alignment.py --parts D --artefacts "baseline drift"
	$(PYTHON) scripts/10_alignment.py --parts D --artefacts "intensity envelope"
	$(PYTHON) scripts/10_alignment.py --parts D --artefacts "additive noise" \
	          --adaptive --suffix " [adaptive]"
	$(PYTHON) scripts/10c_confounders.py
	$(PYTHON) scripts/10b_alignment_figure.py

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf data/processed/* results/* figures/*
