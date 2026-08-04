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


def main() -> int:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    workbook = RAW / "Resultados Informe.xlsx"
    references = RAW / "references.txt"
    missing = [p.name for p in (workbook, references) if not p.exists()]
    if missing:
        print(f"missing from data/raw: {', '.join(missing)}\n", file=sys.stderr)
        print("These files are the property of the laboratory that produced the",
              file=sys.stderr)
        print("characterisation study and are not redistributed here. See the",
              file=sys.stderr)
        print("Data section of the README for how to request them.\n",
              file=sys.stderr)
        print("Steps 2 to 8 do not need them: they read data/processed, which",
              file=sys.stderr)
        print("is included, so the analysis reproduces without this step.",
              file=sys.stderr)
        return 1

    spectra = dataset.load_workbook_spectra(workbook, references)

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
