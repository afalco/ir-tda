# Topological characterisation of shoe-sole materials from IR/ATR spectra

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21792711.svg)](https://doi.org/10.5281/zenodo.21792711)

Persistent homology applied to the infrared spectra of 39 composite materials
(thermoplastic polyurethane, polyurethane, thermoplastic rubber, EVA, PVC and
natural rubber) used in footwear soles.

The method follows the persistence-image workflow developed by Frahi, Falcó,
Chinesta and co-workers for rough surfaces, elastodynamic modes and robot
trajectories, adapted here to one-dimensional vibrational spectra.

Everything needed to reproduce the analysis is in this repository: the
processed spectra and thermogravimetric curves, the material labels, the
library, the scripts and the figures they produce. The raw instrument workbook
is not redistributed — see [Data](#data).

## Results at a glance

See [`docs/REPORT.md`](docs/REPORT.md) for the full discussion.

| Representation | ARI | AMI |
|---|---|---|
| Persistence image (PI) | 0.27 | 0.37 |
| Topological fingerprint image (TFI) | 0.37 | 0.52 |
| Sliced Wasserstein + average linkage | 0.28 | 0.28 |
| Baseline: raw spectra, no topology | 0.56 | 0.66 |

*k*-means with *k* = 3 on the 35 TPU / PUR / TR samples, scored against the
material families, which are never used to build the clusters.

Predicting thermogravimetric behaviour from the spectrum (leave-one-out *Q*²):

| Target | PI | TFI | Family mean | Raw spectra |
|---|---|---|---|---|
| 5 % mass-loss temperature | 0.07 | **0.77** | 0.71 | 0.58 |
| 50 % mass-loss temperature | 0.01 | 0.46 | **0.80** | 0.56 |
| Residue at 800 °C | −0.37 | 0.49 | 0.71 | **0.88** |

Three findings drive the report:

1. **On clean data the raw spectra are hard to beat.** All 39 spectra come from
   one instrument on one calibration, so they already match point by point and
   the invariances that persistence buys are invariances to variation this data
   set does not contain.
2. **Under a wavenumber miscalibration the ordering reverses.** At a ±16 cm⁻¹
   shift the raw spectra identify only 68 % of the materials while the
   persistence image identifies 100 %, because the diagram of a
   one-dimensional filtration does not depend on the parametrisation of the
   axis at all. Topology is the right tool when spectra are pooled across
   instruments, laboratories or calibrations.
3. **On one real task the topological descriptor wins outright.** The
   fingerprint image predicts the onset of thermal degradation with *Q*² = 0.77
   against 0.58 for the raw spectrum, a gap whose paired-bootstrap interval
   excludes zero. The onset is set by the labile minor constituents, which show
   up as moderately prominent bands — exactly what a prominence-weighted
   descriptor is good at. Filler content, encoded in absorbance amplitude, it
   cannot capture at all.

Separating the two polyurethanes, which infrared alone cannot do (ARI):

| Representation | TPU vs PUR |
|---|---|
| Infrared, TFI | 0.10 |
| DTG, TFI | 0.49 |
| Infrared + DTG, TFI | 0.49 |
| Control: raw DTG curve | **0.60** |
| Control: the six conventional TG scalars | 0.22 |

Combining the modalities adds nothing over the thermal one — a negative result,
reported as such. But the persistence image of a DTG curve does separate them
far better than the onset/peak/residue scalars a thermal analyst normally
reports, so the shape of a decomposition profile carries information those
scalars discard.

![robustness](figures/04_robustness.png)

## Layout

```
data/raw/          where the instrument workbook goes; not redistributed
data/processed/    resampled spectra, TG curves and labels — enough to rerun
                   steps 2 to 8 without the workbook
src/irtda/         the library
scripts/           the eight pipeline steps
tests/             unit tests for the persistence and image code
results/           tables produced by the pipeline
figures/           figures produced by the pipeline
docs/REPORT.md     methodology and discussion of the results
environment.yml    conda environment
requirements.txt   pip equivalent
```

## Reproducing the analysis

### Environment

With conda (recommended — it pins the interpreter as well as the libraries):

```bash
git clone https://github.com/afalco/ir-tda.git
cd ir-tda
conda env create -f environment.yml
conda activate ir-tda
```

`environment.yml` installs the package itself in editable mode, so `import
irtda` works from anywhere in the environment. If you prefer mamba, substitute
`mamba env create -f environment.yml`.

With pip and a virtual environment instead:

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

The scripts also run without installing the package at all — each one prepends
`src/` to `sys.path` — so `pip install -r requirements.txt` on its own is
enough to reproduce the figures. Installing is only needed to import `irtda`
from a notebook or from your own code.

Dependencies are the scientific-Python core: numpy, scipy, pandas,
scikit-learn, matplotlib and openpyxl. There is deliberately no topology
library among them, since the persistence computation is implemented directly
(see [Method](#method)).

### Running

```bash
make                      # runs all eight steps, about two minutes
make test                 # 13 unit tests
make clean                # removes everything the pipeline generates
```

or step by step:

```bash
python scripts/01_extract_spectra.py      # workbook  -> spectra + labels
python scripts/02_compute_persistence.py  # spectra   -> diagrams, images, distances
python scripts/03_clustering.py           # features  -> clusters and scores
python scripts/04_robustness.py           # stability under measurement artefacts
python scripts/05_sensitivity.py          # hyper-parameter sweep
python scripts/06_thermal_targets.py      # TG/DTG curves -> processing-window targets
python scripts/07_property_regression.py  # topology -> property prediction
python scripts/08_multimodal.py           # combining the IR and DTG modalities
```

Every step writes plain `.csv` / `.npy` / `.png` files, so intermediate results
can be inspected without rerunning what precedes them.

## Method

Each spectrum is treated as a filtration function on the wavenumber axis. The
degree-0 persistence diagram of its superlevel-set filtration pairs every
absorption band with the level at which it merges into a stronger neighbour, so
that **the persistence of a feature is the topological prominence of a band**:
how far the absorbance has to fall before that band stops being a separate
component. This is exactly the local-maximum/local-minimum pairing described in
*Tape surfaces characterization with persistence images*; for a signal sampled
on a line the Vietoris–Rips construction used there for point clouds collapses
to a union-find computation, which is implemented directly in
`src/irtda/persistence.py` and needs no external topology library.

Diagrams are then vectorised in two ways:

- **PI**, the classical persistence image over `(birth, lifetime)`, weighted by
  a linear ramp in the lifetime and integrated exactly over each pixel;
- **TFI**, a *topological fingerprint image* over `(wavenumber, lifetime)`
  introduced here, which keeps the prominence-based robustness of persistence
  but restores the band positions.

The second one exists because the invariance that makes persistence attractive
for rough surfaces is a liability for spectroscopy: a persistence diagram does
not change if the wavenumber axis is stretched or permuted, yet a band at
1730 cm⁻¹ is an ester carbonyl whatever its height. Replacing the birth
coordinate by the wavenumber of the generating maximum recovers that
information, and it raises the ARI on the TPU/PUR/TR subset from 0.27 to 0.37.

Diagrams are also compared directly with the sliced Wasserstein distance, used
for the dendrogram and the MDS embedding.

## Data

### What is included, and what is not

**Not redistributed.** `data/raw/Resultados Informe.xlsx` and
`data/raw/references.txt` come from the characterisation study
*Caracterización de suelas de calzado*, carried out for a third party by the
Grupo de Investigación de Procesado y Pirólisis de Polímeros, Instituto
Universitario de Ingeniería de Procesos Químicos, Universidad de Alicante
(2021). They are not ours to publish: the workbook holds the full IR and
thermogravimetric traces of 39 commercial compounds and the reference list
names each supplier and trade grade. Requests should go to the authors of that
study.

**Included.** `data/processed/` holds everything the analysis actually needs
downstream: the spectra resampled onto a common grid (`spectra.npz`), the TG
curves (`tg_curves.npz`), the material labels (`labels.csv`) and the thermal
targets (`targets.csv`). Only steps 1 and 6 read the workbook; steps 2 to 8
read `data/processed`, so **the results reproduce in full without it**. Step 1
stops with an explanatory message if the workbook is absent.

For reference, the workbook holds one worksheet per material: columns A and B
are the wavenumber [cm⁻¹] and the ATR absorbance, and a further block is the
thermogravimetric run on the same specimen. The EGA/Py/GC/MS results in the
same worksheets are not used.

`references.txt` maps each reference number to its description and supplier;
the material family is derived from the description by the patterns in
`src/irtda/dataset.py`. The derived families are in `data/processed/labels.csv`,
which is included.

Sheet `H4965` carries no reference number and is mapped to reference 47
(TR "H49 65", Ruiz Alejos) on the strength of the grade name; this is the one
label in the data set assigned by inference rather than read off the list.

## References

1. T. Frahi, C. Argerich, M. Yun, A. Falcó, A. Barasinski, F. Chinesta.
   *Tape surfaces characterization with persistence images*.
   AIMS Materials Science 7(4):364–380, 2020.
2. T. Frahi, F. Chinesta, A. Falcó, A. Badias, E. Cueto, H. Y. Choi, M. Han,
   J.-L. Duval. *Empowering Advanced Driver-Assistance Systems from Topological
   Data Analysis*. Mathematics 9:634, 2021.
3. T. Frahi, A. Falcó, B. Vinh Mau, J. L. Duval, F. Chinesta. *Empowering
   Advanced Parametric Modes Clustering from Topological Data Analysis*.
   Applied Sciences 11:6554, 2021.
4. T. Frahi, A. Sancarlos, M. Galle, X. Beaulieu, A. Chambard, A. Falcó,
   E. Cueto, F. Chinesta. *Monitoring Weeder Robots and Anticipating Their
   Functioning by Using Advanced Topological Data Analysis*.
   Frontiers in Artificial Intelligence 4:761123, 2021.
5. M. Carrière, M. Cuturi, S. Oudot. *Sliced Wasserstein Kernel for Persistence
   Diagrams*. ICML, 2017.

## Citing this work

The release that backs the accompanying manuscript is archived at
<https://doi.org/10.5281/zenodo.21792711>. `CITATION.cff` carries the metadata,
so GitHub's *Cite this repository* button produces a correct entry.

```
Falcó, A. (2026). ir-tda: topological descriptors of infrared spectra for
formulated polymer compounds (v1.0.0). Zenodo.
https://doi.org/10.5281/zenodo.21792711
```

## Licence

The code is MIT, see `LICENSE`. Note that the data in `data/raw/` came from a
characterisation study carried out for a third party and its redistribution is
subject to the terms agreed with its owners; the licence file states this.
