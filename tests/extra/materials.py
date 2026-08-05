"""Tests for practical material registry and file-loading helpers."""

import json

import numpy as np
import pytest

from PyMieSim import (
    MaterialInfo,
    TabulatedMaterial,
    available_materials,
    load_material,
    load_tabulated,
    material_info,
    validate_refractive_indices,
    validate_material,
    validate_wavelength,
    ureg,
)


def test_builtin_registry_exposes_curated_models():
    assert "BK7" in available_materials("sellmeier")
    assert "silver" in available_materials("tabulated")

    info = material_info("BK7")
    assert isinstance(info, MaterialInfo)
    assert info.kind == "sellmeier"
    assert info.wavelength_min < info.wavelength_max


def test_named_material_loader_has_explicit_range_policy():
    material = load_material("BK7")
    assert material.allow_extrapolation is False

    silver = load_material("silver", extrapolation="linear")
    assert silver.allow_extrapolation is True

    with pytest.raises(ValueError, match="outside the model range"):
        validate_wavelength(100 * ureg.nanometer, material)

    validate_material(material)


def test_csv_loader_supports_n_and_k_columns(tmp_path):
    path = tmp_path / "custom.csv"
    path.write_text("wavelength_nm,n,k\n500,1.5,0.01\n600,1.6,0.02\n", encoding="utf-8")

    material = load_tabulated(path)

    assert isinstance(material, TabulatedMaterial)
    np.testing.assert_allclose(material.wavelengths.to("nanometer").magnitude, [500, 600])
    assert material.get_refractive_index(550 * ureg.nanometer) == pytest.approx(1.55 + 0.015j)


def test_json_loader_supports_array_format_and_medium_mode(tmp_path):
    path = tmp_path / "custom.json"
    path.write_text(
        json.dumps({"wavelength": [0.5, 0.6], "n": [1.4, 1.5], "k": [0, 0]}),
        encoding="utf-8",
    )

    medium = load_tabulated(path, medium=True, wavelength_unit="micrometer")

    assert medium.get_refractive_index(550 * ureg.nanometer) == pytest.approx(1.45)


def test_material_validation_rejects_gain_and_bad_wavelengths():
    with pytest.raises(ValueError, match="gain"):
        validate_refractive_indices([1.5 - 0.01j])

    with pytest.raises(ValueError, match="strictly increasing"):
        load_tabulated_data = [500, 500]
        from PyMieSim.materials import validate_tabulated_data

        validate_tabulated_data(load_tabulated_data, [1.5, 1.5])
