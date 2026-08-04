#!/usr/bin/env python3
"""Step 6 -- extract the thermogravimetric targets and the hardness grades.

The workbook records a TG run next to each IR spectrum, on the same specimen.
From it we derive the quantities that bound the thermal processing window of a
compound: the 5, 10 and 50 percent mass-loss temperatures, the temperature of
maximum mass-loss rate, the residue at 800 C and the number of decomposition
steps. Shore hardness grades are parsed from the supplier descriptions.

Writes:
    data/processed/targets.csv
    figures/06_tg_curves.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from irtda import plotting, thermal  # noqa: E402

RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(PROCESSED / "labels.csv")
    xlsx = RAW / "Resultados Informe.xlsx"

    targets = thermal.build_targets(xlsx, labels)
    targets.to_csv(PROCESSED / "targets.csv", index=False)

    available = targets[list(thermal.TG_TARGETS)].notna().all(axis=1).sum()
    print(f"thermal targets available for {available} of {len(targets)} samples")
    print(targets[list(thermal.TG_TARGETS)].describe().round(1).to_string())

    hardness = targets["hardness_raw"].notna().sum()
    print(f"\nhardness grade parsed for {hardness} of {len(targets)} samples")
    print(targets.groupby("hardness_scale")["hardness_raw"]
          .agg(["count", "min", "max"]).to_string())

    # -- figure: TG and DTG curves, coloured by family ----------------------
    curves = thermal.read_tg_curves(xlsx)
    family_of = dict(zip(labels["sheet"].astype(str), labels["family"]))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for sheet, (T, w) in curves.items():
        colour = plotting.FAMILY_COLOURS.get(family_of.get(sheet, "OTHER"), "#7f7f7f")
        grid = np.linspace(max(T.min(), 40.0), min(T.max(), 800.0), 600)
        w_grid = np.interp(grid, T, w)
        axes[0].plot(grid, w_grid, lw=0.8, color=colour, alpha=0.75)
        axes[1].plot(grid, -np.gradient(w_grid, grid), lw=0.8, color=colour, alpha=0.75)

    axes[0].set_xlabel("temperature [$^\\circ$C]")
    axes[0].set_ylabel("residual mass [% w/w]")
    axes[0].set_title("TG")
    axes[1].set_xlabel("temperature [$^\\circ$C]")
    axes[1].set_ylabel("mass-loss rate [% $^\\circ$C$^{-1}$]")
    axes[1].set_title("DTG")

    handles = [plt.Line2D([], [], color=c, lw=2, label=f)
               for f, c in plotting.FAMILY_COLOURS.items()
               if f in set(family_of.values())]
    axes[1].legend(handles=handles, fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)

    fig.suptitle("Thermogravimetric curves of the 39 compounds, by material family")
    fig.tight_layout()
    fig.savefig(FIGURES / "06_tg_curves.png", dpi=180)
    plt.close(fig)

    print(f"\nwritten to {PROCESSED / 'targets.csv'} and "
          f"{FIGURES / '06_tg_curves.png'}")


if __name__ == "__main__":
    main()
