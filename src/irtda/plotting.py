"""Figures for the topological analysis of IR spectra."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from .persistence import lifetime_diagram

__all__ = [
    "plot_spectra_by_family",
    "plot_diagram",
    "plot_persistence_image",
    "plot_dendrogram",
    "plot_embedding",
    "plot_contingency",
    "FAMILY_COLOURS",
]

FAMILY_COLOURS = {
    "TPU": "#1f77b4",
    "PUR": "#d62728",
    "TR": "#2ca02c",
    "EVA": "#ff7f0e",
    "PVC": "#9467bd",
    "RUBBER": "#8c564b",
    "LEATHER": "#e377c2",
    "OTHER": "#7f7f7f",
}


def _colour(family: str) -> str:
    return FAMILY_COLOURS.get(family, "#7f7f7f")


def plot_spectra_by_family(wavenumber, intensity, families, path):
    """One panel per material family, all its normalised spectra overlaid."""
    unique = sorted(set(families))
    n = len(unique)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.6 * nrows),
                             sharex=True, squeeze=False)

    for ax, family in zip(axes.ravel(), unique):
        idx = [i for i, f in enumerate(families) if f == family]
        for i in idx:
            ax.plot(wavenumber, intensity[i], lw=0.6, color=_colour(family), alpha=0.75)
        ax.set_title(f"{family} (n = {len(idx)})", fontsize=10)
        ax.invert_xaxis()
    for ax in axes.ravel()[n:]:
        ax.axis("off")

    fig.supxlabel("wavenumber [cm$^{-1}$]")
    fig.supylabel("normalised absorbance")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_diagram(wavenumber, intensity, diagram, title, path, n_annotate=8):
    """Spectrum, its persistence diagram and its lifetime diagram, side by side."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))

    axes[0].plot(wavenumber, intensity, lw=0.7, color="#333333")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("wavenumber [cm$^{-1}$]")
    axes[0].set_ylabel("normalised absorbance")
    axes[0].set_title("spectrum")

    persistence = diagram[:, 1] - diagram[:, 0]
    lo = min(diagram[:, 0].min(), diagram[:, 1].min())
    hi = max(diagram[:, 0].max(), diagram[:, 1].max())
    axes[1].plot([lo, hi], [lo, hi], color="#bbbbbb", lw=1)
    axes[1].scatter(diagram[:, 0], diagram[:, 1], s=8, c=persistence,
                    cmap="viridis", zorder=3)
    axes[1].set_xlabel("birth")
    axes[1].set_ylabel("death")
    axes[1].set_title(f"persistence diagram ({len(diagram)} features)")

    lt = lifetime_diagram(diagram)
    axes[2].scatter(lt[:, 0], lt[:, 1], s=8, c=lt[:, 1], cmap="viridis")
    axes[2].set_xlabel("birth")
    axes[2].set_ylabel("lifetime")
    axes[2].set_title("lifetime diagram")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_persistence_image(imager, images, labels, path, max_panels=12,
                           title="persistence images"):
    """Grid of persistence images, one per sample."""
    n = min(len(images), max_panels)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 2.4 * nrows),
                             squeeze=False)

    extent = [*imager.birth_range, *imager.lifetime_range]
    for ax, k in zip(axes.ravel(), range(n)):
        img = images[k].reshape(imager.resolution, imager.resolution)
        ax.imshow(img.T, origin="lower", extent=extent, aspect="auto", cmap="viridis")
        ax.set_title(labels[k], fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.ravel()[n:]:
        ax.axis("off")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_dendrogram(distances, labels, families, path, method="average"):
    """Dendrogram of the diagram distance matrix, leaves coloured by family."""
    z = linkage(squareform(distances, checks=False), method=method)

    fig, ax = plt.subplots(figsize=(11, 5))
    # Uniform branch colour: the only colour coding in this figure is the
    # material family of the leaves.
    dendrogram(z, labels=labels, ax=ax, color_threshold=0,
               above_threshold_color="#999999")
    family_of = dict(zip(labels, families))
    for tick in ax.get_xmajorticklabels():
        tick.set_color(_colour(family_of[tick.get_text()]))
        tick.set_fontsize(8)

    ax.set_ylabel(f"{method} linkage on sliced Wasserstein distance")
    handles = [plt.Line2D([], [], color=c, lw=3, label=f)
               for f, c in FAMILY_COLOURS.items() if f in set(families)]
    ax.legend(handles=handles, fontsize=8, ncol=len(handles), loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_embedding(coords, labels, families, path, title):
    """Two-dimensional embedding, annotated with sample identifiers."""
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for family in sorted(set(families)):
        idx = [i for i, f in enumerate(families) if f == family]
        ax.scatter(coords[idx, 0], coords[idx, 1], s=70, color=_colour(family),
                   label=family, edgecolor="white", linewidth=0.8, zorder=3)
    for i, name in enumerate(labels):
        ax.annotate(name, coords[i], fontsize=7, xytext=(4, 4),
                    textcoords="offset points", color="#444444")
    ax.set_title(title)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_contingency(table, path, title):
    """Heat map of a family-versus-cluster contingency table."""
    fig, ax = plt.subplots(figsize=(1.0 * table.shape[1] + 3, 0.55 * table.shape[0] + 2))
    im = ax.imshow(table.values, cmap="Blues", aspect="auto")
    ax.set_xticks(range(table.shape[1]), table.columns, fontsize=9)
    ax.set_yticks(range(table.shape[0]), table.index, fontsize=9)
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            v = table.values[i, j]
            if v:
                ax.text(j, i, str(v), ha="center", va="center",
                        color="white" if v > table.values.max() / 2 else "#333333")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
