"""Correctness checks for the persistence computation and its vectorisation."""

from pathlib import Path

import numpy as np
import pytest

from irtda import clustering, features, images, persistence


def test_single_peak_has_one_finite_class():
    """A unimodal signal has exactly one connected component throughout."""
    x = np.linspace(-3, 3, 201)
    f = np.exp(-(x**2))
    diagram = persistence.superlevel_persistence(f)
    assert len(diagram) == 1
    # The essential class spans the full range of the signal.
    assert persistence.persistence_of(diagram)[0] == pytest.approx(f.max() - f.min())


def test_two_peaks_pair_by_prominence():
    """Two peaks separated by a valley give one finite class of known prominence."""
    #        peak 1      valley     peak 2
    f = np.array([0.0, 1.0, 0.4, 0.8, 0.0])
    diagram = persistence.superlevel_persistence(f)
    lifetimes = np.sort(persistence.persistence_of(diagram))
    # The lower peak (0.8) dies at the valley (0.4): prominence 0.4.
    # The higher peak (1.0) is essential: prominence 1.0.
    assert lifetimes == pytest.approx([0.4, 1.0])


def test_number_of_classes_matches_number_of_local_maxima():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 8 * np.pi, 2000)
    f = np.sin(x) + 0.3 * np.sin(3.7 * x) + 0.01 * rng.normal(size=x.size)
    diagram = persistence.superlevel_persistence(f)
    interior = f[1:-1]
    n_maxima = int(((interior > f[:-2]) & (interior > f[2:])).sum())
    # Every local maximum generates exactly one class (boundary maxima included).
    assert n_maxima <= len(diagram) <= n_maxima + 2


def test_diagram_is_invariant_under_reparametrisation():
    """Persistence of a 1-D signal does not depend on the sampling positions."""
    x = np.linspace(0, 10, 500)
    f = np.sin(x) + 0.5 * np.sin(2.3 * x)
    stretched = np.interp(np.linspace(0, 10, 500), x**1.2 / 10**0.2, f)

    a = persistence.superlevel_persistence(f)
    b = persistence.superlevel_persistence(stretched)
    top_a = np.sort(persistence.persistence_of(a))[-3:]
    top_b = np.sort(persistence.persistence_of(b))[-3:]
    assert top_a == pytest.approx(top_b, abs=0.05)


def test_pruning_keeps_the_most_persistent_features():
    f = np.array([0.0, 1.0, 0.4, 0.8, 0.0, 0.81, 0.79, 0.9, 0.0])
    diagram = persistence.superlevel_persistence(f)
    pruned = persistence.prune_diagram(diagram, min_persistence=0.05)
    kept = persistence.persistence_of(pruned)
    assert (kept > 0.05).all()
    assert len(pruned) < len(diagram)


def test_birth_locations_point_at_the_maxima():
    f = np.array([0.0, 1.0, 0.4, 0.8, 0.0])
    _, born_at, _ = persistence.superlevel_persistence(f, return_locations=True)
    assert set(born_at.tolist()) == {1, 3}


def test_persistence_image_is_translation_covariant():
    """Shifting a diagram along the birth axis shifts its image, not its mass."""
    diagram = np.array([[0.0, 0.5], [0.1, 0.9]])
    shifted = diagram + 0.2

    # Generous padding so that no Gaussian mass falls outside the grid.
    imager = images.PersistenceImager(resolution=16, sigma=0.02, padding=0.5)
    imager.fit([persistence.lifetime_diagram(d) for d in (diagram, shifted)])
    a = imager.transform_one(persistence.lifetime_diagram(diagram))
    b = imager.transform_one(persistence.lifetime_diagram(shifted))
    assert a.sum() == pytest.approx(b.sum(), rel=1e-6)
    assert not np.allclose(a, b)


def test_persistence_image_mass_matches_total_weight():
    """With a narrow kernel well inside the grid, the image integrates the weights."""
    diagram = np.array([[0.0, 1.0], [0.2, 0.6]])
    lt = persistence.lifetime_diagram(diagram)
    imager = images.PersistenceImager(resolution=64, sigma=0.01, padding=0.3)
    imager.fit([lt])
    expected = (lt[:, 1] / lt[:, 1].max()).sum()
    assert imager.transform_one(lt).sum() == pytest.approx(expected, rel=1e-3)


def test_wasserstein_is_zero_for_identical_diagrams():
    diagram = np.array([[0.0, 1.0], [0.2, 0.6], [0.3, 0.35]])
    assert clustering.wasserstein(diagram, diagram) == pytest.approx(0.0, abs=1e-9)
    assert clustering.sliced_wasserstein(diagram, diagram) == pytest.approx(0.0, abs=1e-9)


def test_wasserstein_matches_a_hand_computed_case():
    """One point moved by a known amount, everything else identical."""
    a = np.array([[0.0, 1.0]])
    b = np.array([[0.0, 1.3]])
    assert clustering.wasserstein(a, b, order=2) == pytest.approx(0.3, rel=1e-6)


def test_distance_matrix_is_a_metric_shape():
    rng = np.random.default_rng(1)
    diagrams = [np.sort(rng.random((5, 2)), axis=1) for _ in range(4)]
    d = clustering.distance_matrix(diagrams)
    assert d.shape == (4, 4)
    assert np.allclose(d, d.T)
    assert np.allclose(np.diag(d), 0.0)
    # Triangle inequality on every triple.
    for i in range(4):
        for j in range(4):
            for k in range(4):
                assert d[i, j] <= d[i, k] + d[k, j] + 1e-9


def test_noise_estimate_recovers_the_injected_level():
    rng = np.random.default_rng(2)
    x = np.linspace(0, 20, 5000)
    clean = np.sin(x)
    sigma = 0.02
    noisy = clean + sigma * rng.normal(size=x.size)
    assert features.estimate_noise(noisy) == pytest.approx(sigma, rel=0.1)


def test_multi_config_matches_single_config():
    rng = np.random.default_rng(3)
    wavenumber = np.linspace(500, 4000, 800)
    spectra = np.abs(rng.normal(size=(3, 800)).cumsum(axis=1))
    spectra /= spectra.max(axis=1, keepdims=True)

    a = features.DescriptorConfig(min_persistence=1e-3, adaptive_threshold=False)
    b = features.DescriptorConfig(min_persistence=1e-2, adaptive_threshold=False)
    together = features.compute_diagrams_multi(spectra, wavenumber, [a, b])
    for config, expected in zip((a, b), together):
        got = features.compute_diagrams(spectra, wavenumber, config)
        for x, y in zip(got.diagrams, expected.diagrams):
            assert np.allclose(x, y)


def test_superlevel_diagram_is_expressed_in_spectrum_units():
    """Birth is the band height and death is the merging level, with d < b.

    This pins the sign convention stated in Eqs. (1)-(2) of the paper: the
    superlevel-set diagram lives in the units of ``f`` itself, so its points lie
    strictly below the diagonal and the lifetime is ``ell = b - d > 0``.
    """
    #             peak 1      valley     peak 2
    f = np.array([0.0, 1.0, 0.4, 0.8, 0.0])
    diagram = persistence.superlevel_persistence(f)

    births, deaths = diagram[:, 0], diagram[:, 1]
    assert (deaths < births).all()                      # strictly below diagonal
    assert set(np.round(births, 12)) == {1.0, 0.8}      # births are band heights
    # The lower band (height 0.8) dies at the valley; the essential one at min f.
    assert dict(zip(np.round(births, 12), np.round(deaths, 12))) == {0.8: 0.4, 1.0: 0.0}
    assert persistence.persistence_of(diagram) == pytest.approx([1.0, 0.4])


def test_superlevel_and_sublevel_use_opposite_conventions():
    """``superlevel_persistence(f)`` is the reflection of ``sublevel_persistence(-f)``."""
    rng = np.random.default_rng(7)
    f = np.abs(rng.normal(size=400).cumsum())
    f /= f.max()

    sup = persistence.superlevel_persistence(f)
    sub = persistence.sublevel_persistence(-f)

    assert np.allclose(sup, -sub)
    assert (sup[:, 1] < sup[:, 0]).all()      # superlevel: d < b
    assert (sub[:, 1] > sub[:, 0]).all()      # sublevel:   d > b
    assert np.allclose(persistence.persistence_of(sup), persistence.persistence_of(sub))


def test_lifetime_is_non_negative_for_both_filtrations():
    x = np.linspace(0, 12, 900)
    f = np.sin(x) + 0.4 * np.sin(2.7 * x)

    for diagram in (persistence.superlevel_persistence(f),
                    persistence.sublevel_persistence(f)):
        lt = persistence.lifetime_diagram(diagram)
        assert (lt[:, 1] >= 0).all()
        assert np.allclose(lt[:, 0], diagram[:, 0])


def test_fingerprint_lifetimes_are_positive_prominences():
    wavenumber = np.linspace(500, 4000, 500)
    f = np.zeros(500)
    for centre, height in ((1730.0, 1.0), (2900.0, 0.6), (1100.0, 0.35)):
        f += height * np.exp(-((wavenumber - centre) ** 2) / (2 * 20.0**2))
    f /= f.max()

    diagram, born_at, _ = persistence.superlevel_persistence(f, return_locations=True)
    tfi = persistence.position_lifetime_diagram(diagram, born_at, wavenumber)

    assert (tfi[:, 1] > 0).all()
    assert tfi[:, 0].min() >= wavenumber.min() and tfi[:, 0].max() <= wavenumber.max()
    # The most prominent feature sits at the strongest band.
    assert tfi[np.argmax(tfi[:, 1]), 0] == pytest.approx(1730.0, abs=15.0)


def test_sliced_wasserstein_is_invariant_under_the_sign_convention():
    """Changing convention reflects every diagram, which no distance should see."""
    rng = np.random.default_rng(11)
    f = np.abs(rng.normal(size=600).cumsum()); f /= f.max()
    g = np.abs(rng.normal(size=600).cumsum()); g /= g.max()

    a, b = (persistence.superlevel_persistence(x) for x in (f, g))
    assert clustering.sliced_wasserstein(a, b) == pytest.approx(
        clustering.sliced_wasserstein(-a, -b), rel=1e-9
    )
    assert clustering.wasserstein(a, b) == pytest.approx(
        clustering.wasserstein(-a, -b), rel=1e-9
    )


def test_descriptions_are_translated_into_english():
    """Every Spanish batch description of the study renders in English."""
    from irtda import descriptions

    cases = {
        "5 planchas PUR negro BK sin lavar": "5 black PUR sheets, BK, unwashed",
        "8 pares TR marrón CLO/7922/65": "8 brown TR pairs, CLO/7922/65",
        "4 Planchas de Eva amarillo": "4 yellow EVA sheets",
        "9 muestras PVC negro": "9 black PVC samples",
        "3 planchas caucho marrón LATEX": "3 brown rubber sheets, LATEX",
        "6 planchas EVA negro MII relleno": "6 black EVA sheets, MII, filled",
        # Everything after "marcada como" is a quoted trade name, kept verbatim.
        "2 Planchas C/gris marcada como Hi-react PU Pikolinos":
            "2 grey sheets, labelled Hi-react PU Pikolinos",
    }
    for spanish, english in cases.items():
        assert descriptions.translate(spanish) == english


def test_translation_refuses_an_unknown_spanish_word():
    """A description added later cannot reach a figure untranslated."""
    from irtda import descriptions

    with pytest.raises(descriptions.UnknownTerm):
        descriptions.translate("5 planchas TPU verde XYZ")


def test_every_description_in_the_data_set_is_translatable():
    import pandas as pd
    from irtda import descriptions

    path = Path(__file__).resolve().parents[1] / "data/processed/labels.csv"
    if not path.exists():                       # data not extracted yet
        pytest.skip("data/processed/labels.csv not present")
    for text in pd.read_csv(path)["description"]:
        assert descriptions.translate(text)
