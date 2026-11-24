"""Verify conversion from public PyOptik models to the native solver."""

import numpy as np
import pytest
import PyOptik

from PyMieSim import SellmeierMaterial, SellmeierMedium, TabulatedMaterial, TabulatedMedium, from_pyoptik, ureg


@pytest.mark.parametrize("identifier", [
    "specs/SCHOTT-optical/P-BK7",
    "main/SiO2/Malitson",
    "main/H2O/Daimon-19.0C",
    "main/Ag/Johnson",
    "main/Au/Olmon-ev",
    "main/Al/Rakic",
])
def test_catalog_identity_and_response_are_preserved(identifier):
    source = PyOptik.load_material(name=identifier)
    native = from_pyoptik(source)
    assert source.catalog_id == identifier
    assert native.name == identifier
    bounds = source.wavelength_bound.to("meter").magnitude
    queries = np.linspace(bounds[0], bounds[1], 7)[1:-1] * ureg.meter
    expected = source.compute_refractive_index(wavelength=queries, out_of_range="raise", validation="raise")
    actual = [native.get_refractive_index(wavelength=value) for value in queries]
    assert actual == pytest.approx(expected)
    assert not native.allow_extrapolation


@pytest.mark.parametrize("medium", [False, True])
def test_custom_arrays_select_native_type_and_copy_data(medium):
    source = PyOptik.TabulatedMaterial.from_arrays(
        name="measured", wavelength=[500, 600] * ureg.nanometer, n=[1.5, 1.6],
    )
    native = from_pyoptik(source, medium=medium)
    assert isinstance(native, TabulatedMedium if medium else TabulatedMaterial)
    source.n_values[:] = 2
    assert native.get_refractive_index(wavelength=550 * ureg.nanometer) == pytest.approx(1.55)


@pytest.mark.parametrize("medium", [False, True])
def test_formula_conversion_selects_native_type(medium):
    source = PyOptik.load_material(name="main/H2O/Daimon-19.0C")
    native = from_pyoptik(source, medium=medium)
    assert isinstance(native, SellmeierMedium if medium else SellmeierMaterial)


def test_csv_loading_is_owned_by_pyoptik(tmp_path):
    path = tmp_path / "measurement.csv"
    path.write_text("wavelength,n,k\n500,1.5,\n600,1.6,0.02\n")
    source = PyOptik.TabulatedMaterial.from_csv(path=path, wavelength_unit=ureg.nanometer, missing_k=0.0)
    native = from_pyoptik(source)
    assert native.get_refractive_index(wavelength=550 * ureg.nanometer) == pytest.approx(1.55 + 0.01j)


def test_separate_component_grids_preserve_both_linear_curves(tmp_path):
    path = tmp_path / "split.yml"
    path.write_text("DATA:\n- type: tabulated n\n  data: |\n    0.4 1.4\n    0.5 1.8\n    0.7 1.9\n- type: tabulated k\n  data: |\n    0.45 0.1\n    0.55 0.4\n    0.6 0.2\n")
    source = PyOptik.TabulatedMaterial(filename="split", file_path=path)
    native = from_pyoptik(source)
    for wavelength in np.linspace(460, 590, 15) * ureg.nanometer:
        assert native.get_refractive_index(wavelength=wavelength) == pytest.approx(
            source.compute_refractive_index(wavelength=wavelength, out_of_range="raise")
        )


@pytest.mark.parametrize("index", [-1.5, 1.5 - 0.01j])
def test_strict_validation_is_delegated_to_pyoptik(index):
    source = PyOptik.TabulatedMaterial.from_arrays(
        name="sample", wavelength=[500, 600] * ureg.nanometer,
        n=[complex(index).real] * 2, k=[complex(index).imag] * 2,
    )
    with pytest.raises(ValueError, match="positive|gain"):
        from_pyoptik(source)
    native = from_pyoptik(source, validate=False)
    assert native.get_refractive_index(wavelength=550 * ureg.nanometer) == pytest.approx(index)


@pytest.mark.parametrize("validate", [True, False])
def test_medium_extinction_cannot_be_silently_discarded(validate):
    source = PyOptik.TabulatedMaterial.from_arrays(
        name="absorbing", wavelength=[500, 600] * ureg.nanometer, n=[1.5, 1.6], k=[0.01, 0.02],
    )
    with pytest.raises(ValueError, match="lossless medium"):
        from_pyoptik(source, medium=True, validate=validate)


@pytest.mark.parametrize("formula", [False, True])
def test_range_policy_is_applied_to_native_model(formula):
    source = (PyOptik.SellmeierMaterial.from_coefficients(
        name="formula", coefficients=[1.5, 0.1, 1], formula_type=5,
        wavelength_range=[500, 600] * ureg.nanometer,
    ) if formula else PyOptik.TabulatedMaterial.from_arrays(
        name="table", wavelength=[500, 600] * ureg.nanometer, n=[1.5, 1.6],
    ))
    strict = from_pyoptik(source)
    with pytest.raises(IndexError, match="outside|range"):
        strict.get_refractive_index(wavelength=650 * ureg.nanometer)
    permissive = from_pyoptik(source, extrapolation="linear")
    assert permissive.allow_extrapolation
    assert np.isfinite(permissive.get_refractive_index(wavelength=650 * ureg.nanometer))


def test_cubic_interpolation_requires_explicit_resampling():
    source = PyOptik.TabulatedMaterial.from_arrays(
        name="cubic", wavelength=[500, 600] * ureg.nanometer, n=[1.5, 1.6], interpolation="pchip",
    )
    with pytest.raises(ValueError, match="Only linear"):
        from_pyoptik(source)


def test_unsupported_formula_is_rejected():
    source = PyOptik.SellmeierMaterial.from_coefficients(
        name="polynomial", coefficients=[2, 1, 2], formula_type=3,
        wavelength_range=[500, 600] * ureg.nanometer,
    )
    with pytest.raises(ValueError, match="Formula 3 is not supported"):
        from_pyoptik(source)


def test_unsupported_model_and_range_policy_are_rejected():
    with pytest.raises(TypeError, match="Expected a PyOptik"):
        from_pyoptik(object())
    source = PyOptik.TabulatedMaterial.from_arrays(name="sample", wavelength=[0.5, 0.6], n=[1.5, 1.6])
    with pytest.raises(ValueError, match="extrapolation"):
        from_pyoptik(source, extrapolation="clip")
