#!/usr/bin/env python
"""Regression tests for the high-level experiment API."""

import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PyMieSim import LabeledArray, Measure
from PyMieSim.experiment import Setup
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.units import ureg


def _setup():
    source = GaussianSet(
        wavelength=[600, 700] * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=[100, 200] * ureg.nanometer,
        material=[1.4],
        medium=[1.0],
    )
    return Setup(scatterer_set=scatterer, source_set=source)


def test_get_rejects_empty_measure_list():
    with pytest.raises(ValueError, match="At least one measure"):
        _setup().get()


def test_get_rejects_unknown_measure():
    with pytest.raises(ValueError, match="Unknown measure"):
        _setup().get("not_a_measure")


def test_labeled_array_as_numpy_preserves_measure_order_and_shape():
    experiment = _setup()

    values = experiment.get("Qext", "Qsca").as_numpy()

    assert values.shape == (2, 2, 2)
    np.testing.assert_allclose(values[0], experiment.get("Qext").as_numpy())
    np.testing.assert_allclose(values[1], experiment.get("Qsca").as_numpy())


def test_labeled_array_as_numpy_keeps_simulation_shape():
    values = _setup().get("Qsca").as_numpy()

    assert values.shape == (2, 2)


def test_labeled_array_converts_to_numpy_and_dataframe():
    result = _setup().get("Qsca")

    np.testing.assert_allclose(result.as_numpy(), result.to_numpy())
    dataframe = result.as_dataframe()

    assert type(dataframe).__name__ == "DataFrame"
    assert list(dataframe.columns) == [
        "source:wavelength",
        "scatterer:diameter",
        "Qsca",
    ]
    assert dataframe.attrs["measures"] == ("Qsca",)


def test_labeled_array_converts_multiple_measures_to_tabular_columns():
    result = _setup().get("Qext", "Qsca")

    values = result.as_numpy()
    dataframe = result.as_dataframe()

    assert values.shape == (2, 2, 2)
    assert type(dataframe).__name__ == "DataFrame"
    assert list(dataframe.columns) == [
        "source:wavelength",
        "scatterer:diameter",
        "Qext",
        "Qsca",
    ]
    assert len(dataframe) == 4


def test_get_labeled_array_preserves_grid_coordinates():
    experiment = _setup()

    result = experiment.get("Qsca")

    assert isinstance(result, LabeledArray)
    assert result.dims == ("source:wavelength", "scatterer:diameter")
    assert result.shape == (2, 2)
    np.testing.assert_allclose(result.to_numpy(), experiment.get("Qsca").as_numpy())
    np.testing.assert_allclose(result.coords["source:wavelength"], [600e-9, 700e-9])
    assert result.attrs["coordinate_units"]["source:wavelength"] == ureg.meter


def test_labeled_array_repr_is_informative():
    result = _setup().get("Qsca")

    representation = repr(result)

    assert representation.startswith("LabeledArray(")
    assert "shape=(2, 2)" in representation
    assert "source:wavelength" in representation


def test_get_labeled_array_stacks_measures_on_a_named_dimension():
    result = _setup().get("Qext", "Qsca")

    assert result.dims[0] == "measure"
    assert result.shape == (2, 2, 2)
    assert list(result.coords["measure"]) == ["Qext", "Qsca"]


def test_labeled_array_supports_indexing_and_reduction():
    result = _setup().get("Qsca")

    selected = result.isel({"source:wavelength": 0})
    reduced = result.mean("source:wavelength")
    reduced_with_x = result.mean(x="source:wavelength")

    assert selected.dims == ("scatterer:diameter",)
    assert selected.shape == (2,)
    assert reduced.dims == ("scatterer:diameter",)
    np.testing.assert_allclose(reduced.to_numpy(), result.to_numpy().mean(axis=0))
    np.testing.assert_allclose(reduced_with_x.to_numpy(), reduced.to_numpy())


def test_labeled_array_plot_supports_parameter_x_and_measure_y():
    result = _setup().get("Qsca")
    figure, axis = plt.subplots()

    artists = result.plot(
        x="scatterer:diameter",
        y="Qsca",
        ax=axis,
    )

    assert len(artists) == 2
    assert all(artist in axis.lines for artist in artists)
    assert axis.get_legend() is not None
    assert len(axis.get_legend().texts) == 2
    assert all("Qsca |" in text.get_text() for text in axis.get_legend().texts)
    plt.close(figure)


def test_labeled_array_plot_uses_parameter_x_and_measure_y():
    result = _setup().get("Qsca").mean(x="source:wavelength")
    figure, axis = plt.subplots()

    artists = result.plot(
        x="scatterer:diameter",
        y="Qsca",
        ax=axis,
    )

    assert len(artists) == 1
    assert artists[0] in axis.lines
    plt.close(figure)


def test_labeled_array_plot_supports_standard_deviation_band():
    result = _setup().get("Qsca")
    figure, axis = plt.subplots()

    artists = result.plot(
        x="source:wavelength",
        y="Qsca",
        std="scatterer:diameter",
        ax=axis,
    )

    assert len(artists) == 2
    assert artists[0] in axis.lines
    assert artists[1] in axis.collections
    assert axis.get_legend() is not None
    assert axis.get_xlabel() == "source:wavelength"
    assert axis.get_ylabel() == "Qsca"
    plt.close(figure)


def test_labeled_array_plot_supports_remaining_parameter_sweeps():
    values = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
    result = LabeledArray(
        values,
        ["detector:NA", "detector:phi_offset", "detector:cache_numerical_aperture"],
        {
            "detector:NA": [0.1, 0.2],
            "detector:phi_offset": [0, 1, 2],
            "detector:cache_numerical_aperture": [0.05, 0.08, 0.1, 0.12],
        },
        {"measures": ("coupling",)},
        "coupling",
    )
    figure, axis = plt.subplots()

    artists = result.plot(
        x="detector:phi_offset",
        y="coupling",
        std="detector:NA",
        ax=axis,
    )

    assert len(axis.lines) == 4
    assert len(axis.collections) == 4
    assert len(artists) == 8
    assert len(axis.get_legend().texts) == 4
    plt.close(figure)


def test_experiment_exposes_available_measures_and_accepts_enum():
    experiment = _setup()

    assert "Qsca" in experiment.available_measures
    result = experiment.get(Measure.QSCA)

    assert isinstance(result, LabeledArray)
    assert result.name == "Qsca"
    assert "scatterer:diameter" in result.dims
