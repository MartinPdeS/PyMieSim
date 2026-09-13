"""Single and sweep APIs share validation, physical results, and typed metadata."""
import numpy as np
import pytest

from PyMieSim import (
    Experiment, Gaussian, GaussianSet, Measure, PolarizationState, PolarizationSet,
    Simulation, SimulationResult, SimulationResults, Sphere, SphereSet, ureg,
)


def configurations():
    simulation = Simulation(
        scatterer=Sphere(diameter=100 * ureg.nanometer, material=1.4, medium=1.0),
        source=Gaussian(wavelength=600 * ureg.nanometer,
                        polarization=PolarizationState(angle=0 * ureg.degree),
                        optical_power=1e-3 * ureg.watt, numerical_aperture=0.2),
    )
    experiment = Experiment(
        scatterer_set=SphereSet(diameter=[100, 200] * ureg.nanometer, material=[1.4], medium=[1.0]),
        source_set=GaussianSet(wavelength=[600, 700] * ureg.nanometer,
                              polarization=PolarizationSet(angles=0 * ureg.degree),
                              optical_power=[1e-3] * ureg.watt, numerical_aperture=[0.2]),
    )
    return simulation, experiment


@pytest.mark.parametrize('names, message', [
    ((), 'At least one'), (('invalid',), 'Unknown measure'),
    (('coupling',), 'Unknown measure'), (('Qsca', Measure.QSCA), 'only once'),
    (([],), 'Unknown measure'),
])
def test_shared_validation(names, message):
    for model in configurations():
        for method in (model.run, model.get):
            with pytest.raises(ValueError, match=message):
                method(*names)


def test_multiple_single_measures_and_sweep_agree():
    simulation, experiment = configurations()
    raw = simulation.run(Measure.QSCA, 'Csca')
    assert tuple(raw) == ('Qsca', 'Csca')
    single = simulation.get(Measure.QSCA, 'Csca', as_result=True)
    sweep = experiment.run(Measure.QSCA, 'Csca', as_result=True)
    assert isinstance(single, SimulationResults)
    assert isinstance(sweep, SimulationResults)
    for name in raw:
        assert single[name].units == sweep[name].units
        np.testing.assert_allclose(raw[name].magnitude, sweep[name].magnitude[0, 0])
        np.testing.assert_allclose(single[name].magnitude, raw[name].magnitude)


def test_typed_sweep_preserves_coordinates_and_unit_conversion():
    _, experiment = configurations()
    result = experiment.get('Csca', as_result=True)
    assert isinstance(result, SimulationResult)
    converted = result.to('nanometer ** 2')
    assert converted.dims == ('source:wavelength', 'scatterer:diameter')
    np.testing.assert_allclose(converted.coords['source:wavelength'], [600e-9, 700e-9])
    assert converted.coordinate_units['source:wavelength'] == ureg.meter
    labeled = converted.as_labeled_array()
    assert labeled.attrs['units']['Csca'] == ureg.nanometer ** 2
    np.testing.assert_allclose(labeled.as_numpy(), experiment.get('Csca').as_numpy() * 1e18)
    np.testing.assert_allclose(result.as_labeled_array().as_numpy(), experiment.run('Csca').as_numpy())


def test_typed_results_keep_singleton_axes_when_requested():
    _, experiment = configurations()
    result = experiment.run('Qsca', drop_unique_level=False, as_result=True)
    expected = experiment.get('Qsca', drop_unique_level=False)
    assert result.dims == expected.dims
    assert np.shape(result.magnitude) == expected.shape
    assert len(result.dims) > 2


def test_single_result_exports_scalar_labeled_array():
    simulation, _ = configurations()
    result = simulation.run('Csca', as_result=True).to('nanometer ** 2')
    labeled = result.as_labeled_array()
    assert labeled.shape == ()
    assert labeled.dims == ()
    assert labeled.name == 'Csca'
    np.testing.assert_allclose(labeled.as_numpy(), result.magnitude)


@pytest.mark.parametrize('method, field_count', [('get_farfields', 2), ('get_stokes', 4)])
def test_field_overloads_match_native_return_shapes(method, field_count):
    simulation, _ = configurations()
    structured = getattr(simulation, method)(sampling=4, distance=1 * ureg.meter)
    assert len(structured) == field_count + 1
    for field in structured[:-1]:
        assert field.magnitude.shape == (4, 4)
    sampled = getattr(simulation, method)(
        phi=np.array([0., 30.]) * ureg.degree,
        theta=np.array([0., 45.]) * ureg.degree,
        distance=1 * ureg.meter,
    )
    assert len(sampled) == field_count
    for field in sampled:
        assert field.magnitude.shape == (2,)
