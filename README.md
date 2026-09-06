# Topological characterisation of shoe-sole materials from IR/ATR spectra

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21795105.svg)](https://doi.org/10.5281/zenodo.21795105)

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

Six findings drive the report. The last two were established by the controls in
steps 9 to 12 and correct earlier conclusions of this work; they are stated here
rather than in a footnote because they narrow what the descriptor can be claimed
to do.

1. **On clean data the raw spectra are hard to beat.** All 39 spectra come from
   one instrument on one calibration, so they already match point by point and
   the invariances that persistence buys are invariances to variation this data
   set does not contain.
2. **On one task the topological descriptor leads.** The fingerprint image
   predicts the onset of thermal degradation with *Q*² = 0.77 against 0.58 for
   the raw spectrum, 0.67 for a second-derivative chain and 0.67 for a baseline
   correction; the paired-bootstrap interval against the raw spectrum excludes
   zero, and a permutation test puts it far beyond its null (*p* = 0.003). The
   onset is set by the labile minor constituents, which show up as moderately
   prominent bands — exactly what a prominence-weighted descriptor reads. Filler
   content, encoded in absorbance amplitude, it cannot capture at all.
3. **Persistence and a peak table carry the same quantity.** Degree-0
   persistence of a superlevel-set filtration is the topographic prominence a
   peak picker reports: at matched thresholds 98.9 % of the topological features
   coincide with a detected peak, and on 96 % of interior bands the two values
   agree to machine precision. A conventional peak table carrying position and
   prominence reaches *Q*² = 0.72 on the onset, statistically indistinguishable
   from the fingerprint image (Δ*Q*² = 0.05, [−0.05, 0.12]). What persistence
   supplies is that quantity computed exactly and with no detection threshold,
   not information the peak table lacks. Among the conventional attributes it is
   prominence that carries the signal: position alone, intensity, width and area
   are all clearly worse.
4. **The robustness is not topological.** Under a ±16 cm⁻¹ miscalibration the
   untreated raw spectra retrieve 0.79 and the band-based descriptors 0.97, but
   aligning the query by cross-correlation restores the raw spectra to 0.99 and
   a shift-invariant Fourier magnitude reaches 1.00. Across all four artefacts
   the worst case is 0.67 for a baseline correction with SNV and 0.18 for the
   fingerprint image. No general claim of robustness survives.
5. **The Vietoris–Rips pipeline of the literature loses on these data.**
   Implemented as described — each spectral sample a point of the plane, a
   Vietoris–Rips filtration, Betti curves — it reaches *Q*² = 0.14 on the onset
   and ARI 0.24 on the families, against 0.77 and 0.37 for the fingerprint
   image. Its *H*₁ classes predict the onset not at all, and their number varies
   by a factor of twelve with the ratio between the two axes, which nothing in
   the measurement fixes. It also costs 108× the lower-star sweep and has to
   subsample the spectrum to run. The objection it rests on — that a lower-star
   filtration yields no *H*₁ for a 1-D signal — is correct and, here, costs
   nothing.
6. **The two polyurethanes are separable by infrared after all.** A
   Savitzky–Golay second derivative with SNV recovers TPU, PUR and TR exactly by
   *k*-means (ARI 1.000), stable over ten of twelve filter settings and
   attributable neither to the supplier who delivered each batch nor to its
   colour — each family stays in one cluster across two or three different
   suppliers. An earlier version of this work concluded the opposite.

Separating the two polyurethanes (ARI):

| Representation | TPU vs PUR |
|---|---|
| Infrared, TFI | 0.10 |
| DTG, TFI | 0.49 |
| Infrared + DTG, TFI | 0.49 |
| Control: raw DTG curve | 0.60 |
| Control: the six conventional TG scalars | 0.22 |
| **Control: infrared, Savitzky–Golay 2nd derivative** | **1.00** |

Combining the infrared and thermal modalities adds nothing over the thermal one
— a negative result, reported as such. The persistence image of a DTG curve does
separate the two polyurethanes far better than the onset/peak/residue scalars a
thermal analyst normally reports, so the shape of a decomposition profile
carries information those scalars discard. But the cleanest separation of all
comes from the infrared spectrum with an ordinary derivative filter, which is
what the last row records.

The onset prediction expressed as a processing decision rather than as a *Q*²,
with a margin that carries a finite-sample guarantee (jackknife+; realised risk
2.6 % against a 5 % nominal):

| Representation | mean abs. error | certified margin |
|---|---|---|
| Family mean | 10.8 °C | 26.5 °C |
| Raw spectra (SNV) | 11.4 °C | 46.3 °C |
| Savitzky–Golay 2nd derivative | 10.1 °C | 68.9 °C |
| Peaks, position + prominence | **8.9 °C** | **24.3 °C** |
| Fingerprint image | 9.0 °C | 24.6 °C |

Two degrees of moulding temperature returned over knowing the polymer class
(1.96 °C, [1.08, 2.86]) — certified, where an uncalibrated empirical quantile
suggested eight and delivered a 7.7 % risk when asked for 5 %. Note that the
second-derivative chain has the second best typical error and by far the worst
certified margin: the margin is set by the lower tail of the residual, and
identification and prediction are not won by the same representation.

![robustness](figures/04_robustness.png)

## Layout

```
data/raw/          where the instrument workbook goes; not redistributed
data/processed/    resampled spectra, TG curves and labels — enough to rerun
                   steps 2 to 13 without the workbook
src/irtda/         the library: persistence, images, clustering, dataset,
                   features, thermal, plotting, plus descriptions (English
                   rendering of the Spanish batch labels), peaks (conventional
                   peak descriptors, the control of step 9) and preprocess
                   (derivative, baseline, scatter-correction and alignment
                   baselines, the controls of step 10) and rips (the
                   Vietoris-Rips pipeline of the literature, the control of
                   step 13)
scripts/           the sixteen pipeline steps
tests/             unit tests for the persistence, image, description, peak
                   and pre-processing code
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
make                      # runs every step, a few minutes
make test                 # 31 unit tests
make clean                # removes everything the pipeline generates
```

or by stage:

```bash
make extract persistence            # workbook -> spectra, diagrams, images
make cluster robustness sensitivity
make thermal regression multimodal
make peaks                          # controlled comparison against peak tables
make alignment                      # pre-processing and alignment baselines
make window conformal               # the moulding decision and its guarantee
make rips                           # the Vietoris-Rips pipeline of the literature
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
python scripts/09_peak_features.py        # is the fingerprint more than a peak table?
python scripts/09b_peak_figure.py         #   ... and its figure
python scripts/10_alignment.py            # pre-processing and alignment baselines
python scripts/10c_confounders.py         # batch effects behind the derivative result
python scripts/10b_alignment_figure.py    #   ... and its figure
python scripts/11_moulding_window.py      # the prediction as a moulding decision
python scripts/12_conformal_margin.py     # a safety margin with a guarantee
python scripts/13_rips_comparison.py      # the Vietoris-Rips pipeline, head to head
```

Steps 9 and 10 are the long ones. Both take `--parts`, and step 9 also
`--targets` while step 10 takes `--artefacts`, so an expensive stage can be run
one target or one artefact at a time; each invocation adds its rows to the same
table rather than overwriting it.

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
<https://doi.org/10.5281/zenodo.21795105>. `CITATION.cff` carries the metadata,
so GitHub's *Cite this repository* button produces a correct entry, and
`.zenodo.json` supplies the archival metadata so that every future release is
described correctly without editing anything by hand.

```
Falcó, A. (2026). ir-tda: topological descriptors of infrared spectra for
formulated polymer compounds (v1.1.0). Zenodo.
https://doi.org/10.5281/zenodo.21795105
```


## Licence

The code is MIT, see `LICENSE`. Note that the data in `data/raw/` came from a
characterisation study carried out for a third party and its redistribution is
subject to the terms agreed with its owners; the licence file states this.
