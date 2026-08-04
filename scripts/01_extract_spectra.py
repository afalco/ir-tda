#!/usr/bin/env python3
"""Step 1 -- extract the IR/ATR spectra and the material labels from the raw data.

Reads ``data/raw/Resultados Informe.xlsx`` (one worksheet per material
reference) together with ``data/raw/references.txt``, resamples every spectrum
onto a common wavenumber grid and writes:

    data/processed/spectra.npz   wavenumber grid and raw absorbance matrix
    data/processed/labels.csv    sample identifier, description, supplier, family
    figures/01_spectra_by_family.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from irtda import dataset, plotting  # noqa: E402

RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    spectra = dataset.load_workbook_spectra(
        RAW / "Resultados Informe.xlsx", RAW / "references.txt"
    )

    np.savez_compressed(
        PROCESSED / "spectra.npz",
        wavenumber=spectra.wavenumber,
        intensity=spectra.intensity,
        sheet=spectra.metadata["sheet"].to_numpy(),
    )
    spectra.metadata.to_csv(PROCESSED / "labels.csv", index=False)

    normalised = dataset.normalise(spectra.intensity, "minmax")
    plotting.plot_spectra_by_family(
        spectra.wavenumber,
        normalised,
        spectra.metadata["family"].tolist(),
        FIGURES / "01_spectra_by_family.png",
    )

    print(f"{len(spectra)} spectra x {spectra.intensity.shape[1]} points "
          f"({spectra.wavenumber[0]:.0f}-{spectra.wavenumber[-1]:.0f} cm-1)")
    print(spectra.metadata["family"].value_counts().to_string())
    print(f"\nwritten to {PROCESSED}")


if __name__ == "__main__":
    main()
