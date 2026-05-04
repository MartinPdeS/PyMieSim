#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from PyMieSim.units import ureg
from PyMieSim.gui.parsing import parse_material_values, parse_numeric_expression, parse_quantity_expression
from PyMieSim.gui.services import available_measures, build_detector_set, run_experiment


def test_parse_quantity_expression_supports_ranges():
    quantity = parse_quantity_expression("400:800:5", ureg.nanometer)

    assert len(quantity) == 5
    assert quantity.units == ureg.nanometer
    assert np.isclose(quantity.magnitude[-1], 800)


def test_parse_material_values_supports_named_materials():
    material = parse_material_values("silver")

    assert material.__class__.__name__ == "TabulatedMaterial"


def test_parse_numeric_expression_supports_integer_scalars():
    sampling = parse_numeric_expression("200", integer=True)

    assert sampling == 200


def test_available_measures_excludes_coupling_without_detector():
    measures = available_measures("SphereSet", "None")

    assert "coupling" not in measures
    assert "Qsca" in measures


def test_build_detector_set_returns_none_for_detectorless_runs():
    detector = build_detector_set("None", {})

    assert detector is None


def test_run_experiment_returns_serialized_dataframe():
    result = run_experiment(
        source_type="GaussianSet",
        source_values={
            "wavelength": "650",
            "polarization": "0",
            "optical_power": "1e-3",
            "numerical_aperture": "0.2",
        },
        scatterer_type="SphereSet",
        scatterer_values={
            "diameter": "500:700:3",
            "material": "1.4",
            "medium": "1.0",
        },
        detector_type="PhotodiodeSet",
        detector_values={
            "numerical_aperture": "0.2",
            "gamma_offset": "0",
            "phi_offset": "0",
            "sampling": "50",
        },
        measure="Qsca",
    )

    assert result["measure"] == "Qsca"
    assert result["row_count"] == 3
    assert len(result["rows"]) == 3
    assert result["parameter_columns"]