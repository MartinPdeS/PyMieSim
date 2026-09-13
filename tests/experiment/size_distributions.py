"""Physical weighting, distribution moments, and independent-scattering averages."""
import numpy as np
import pytest
from pint import DimensionalityError

from PyMieSim import (
    Experiment, ParticleSizeDistribution, PlaneWaveSet, PolarizationSet,
    SphereSet, InfiniteCylinderSet, SimulationResult, SimulationResults, ureg,
)


def experiment_for(distribution, material=1.5):
    return Experiment(
        source_set=PlaneWaveSet(wavelength=[500, 700] * ureg.nanometer,
                                polarization=PolarizationSet(angles=0 * ureg.degree),
                                amplitude=[1] * ureg.volt / ureg.meter),
        scatterer_set=SphereSet(diameter=distribution.diameters, material=[material], medium=[1.0]),
    )


def average(experiment, distribution, *names, **kwargs):
    return experiment.average_size_distribution(
        distribution, *names, **kwargs,
    )


@pytest.mark.parametrize('weights', [[1, 3], [1e308, 1e308]])
def test_number_fractions_are_normalized_and_inputs_copied(weights):
    sizes = np.array([100., 300.]) * ureg.nanometer
    distribution = ParticleSizeDistribution(sizes, weights)
    assert np.sum(distribution.number_fractions) == pytest.approx(1)
    sizes[0] = 800 * ureg.nanometer
    fractions = distribution.number_fractions
    fractions[:] = 0
    assert distribution.diameters.to('nanometer').magnitude[0] == pytest.approx(100)
    assert distribution.number_fractions.sum() == pytest.approx(1)


@pytest.mark.parametrize('diameters, weights', [
    ([], []), ([0], [1]), ([-1], [1]), ([np.inf], [1]), ([np.nan], [1]),
    ([1, 2], [1]), ([[1, 2]], [[1, 1]]), ([1], [-1]), ([1], [0]),
    ([1], [np.nan]), ([1], [np.inf]), ([1 + 1j], [1]), ([1], [1j]),
])
def test_invalid_discrete_distributions(diameters, weights):
    with pytest.raises(ValueError):
        ParticleSizeDistribution(np.array(diameters) * ureg.nanometer, weights)


def test_length_units_are_required():
    with pytest.raises(TypeError, match='length units'):
        ParticleSizeDistribution([100, 200], [1, 1])
    with pytest.raises(DimensionalityError):
        ParticleSizeDistribution([100, 200] * ureg.second, [1, 1])


@pytest.mark.parametrize('order', [1, 2, 3, 6])
def test_lognormal_quadrature_matches_analytic_moments(order):
    median = 100 * ureg.nanometer
    geometric_std = 1.4
    distribution = ParticleSizeDistribution.lognormal(median, geometric_std, sampling=32)
    actual = np.sum(distribution.number_fractions * distribution.diameters.to('nanometer').magnitude ** order)
    expected = 100 ** order * np.exp(0.5 * order ** 2 * np.log(geometric_std) ** 2)
    assert actual == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize('kwargs', [
    {'geometric_std': 0.9}, {'geometric_std': np.inf}, {'geometric_std': np.nan},
    {'sampling': 0}, {'sampling': 2.5}, {'sampling': True},
    {'median_diameter': -1 * ureg.nanometer}, {'median_diameter': [1, 2] * ureg.nanometer},
])
def test_invalid_lognormal_parameters(kwargs):
    parameters = {'median_diameter': 100 * ureg.nanometer, 'geometric_std': 1.2, **kwargs}
    with pytest.raises(ValueError):
        ParticleSizeDistribution.lognormal(**parameters)


def test_monodisperse_limit_matches_native_sphere():
    distribution = ParticleSizeDistribution.lognormal(100 * ureg.nanometer, 1)
    assert len(distribution.number_fractions) == 1
    experiment = experiment_for(distribution)
    for name in ('Csca', 'Cext', 'Qsca', 'Qext', 'g'):
        result = average(experiment, distribution, name)
        np.testing.assert_allclose(result.as_numpy(), experiment.get(name).as_numpy(), rtol=1e-12)


def test_cross_sections_and_asymmetry_use_distinct_physical_weights():
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 3])
    experiment = experiment_for(distribution)
    csca = experiment.get('Csca').as_numpy()
    asymmetry = experiment.get('g').as_numpy()
    weights = distribution.number_fractions
    results = average(experiment, distribution, 'Csca', 'g', 'Qsca', as_result=True)
    assert isinstance(results, SimulationResults)
    expected_csca = csca @ weights
    expected_g = (csca * asymmetry) @ weights / expected_csca
    expected_area = np.sum(weights * np.pi / 4 * distribution.diameters.magnitude ** 2)
    np.testing.assert_allclose(results['Csca'].magnitude, expected_csca)
    np.testing.assert_allclose(results['g'].magnitude, expected_g)
    np.testing.assert_allclose(results['Qsca'].magnitude, expected_csca / expected_area)
    assert not np.allclose(expected_g, asymmetry @ weights, rtol=1e-3)
    assert not np.allclose(results['Qsca'].magnitude, experiment.get('Qsca').as_numpy() @ weights, rtol=1e-3)


def test_energy_balance_and_metadata_survive_conversion_and_export():
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 3])
    experiment = experiment_for(distribution, material=1.5 + 0.1j)
    result = average(experiment, distribution, 'Cext', 'Csca', 'Cabs')
    np.testing.assert_allclose(result.as_numpy()[0], result.as_numpy()[1] + result.as_numpy()[2], rtol=1e-12)
    assert result.dims == ('measure', 'source:wavelength')
    assert result.attrs['approximation'] == 'independent_scattering'
    assert result.as_dataframe().attrs['approximation'] == 'independent_scattering'
    typed = average(experiment, distribution, 'Csca', as_result=True).to('nanometer ** 2')
    assert isinstance(typed, SimulationResult)
    assert typed.metadata['approximation'] == 'independent_scattering'
    assert typed.as_labeled_array().attrs['assumptions'] == result.attrs['assumptions']
    assert typed.coordinate_units == {'source:wavelength': ureg.meter}
    np.testing.assert_allclose(typed.coords['source:wavelength'], [500e-9, 700e-9])


def test_singleton_axes_and_zero_probability_bin():
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 0])
    experiment = experiment_for(distribution)
    result = average(experiment, distribution, 'Csca', drop_unique_level=False)
    raw = experiment.get('Csca', drop_unique_level=False)
    axis = raw.dims.index('scatterer:diameter')
    np.testing.assert_allclose(result.as_numpy(), np.take(raw.as_numpy(), 0, axis=axis))
    assert 'scatterer:diameter' not in result.dims
    assert len(result.dims) == len(raw.dims) - 1


def test_averaging_needs_no_approximation_argument():
    distribution = ParticleSizeDistribution([100] * ureg.nanometer, [1])
    experiment = experiment_for(distribution)
    result = experiment.average_size_distribution(distribution, 'Csca')
    np.testing.assert_allclose(result.as_numpy(), experiment.get('Csca').as_numpy())


@pytest.mark.parametrize('names', [(), ('a1',), ('Qratio',), ('Cratio',), ('coupling',), ('Csca', 'Csca')])
def test_unsupported_or_ambiguous_measures_are_rejected(names):
    distribution = ParticleSizeDistribution([100] * ureg.nanometer, [1])
    with pytest.raises(ValueError):
        average(experiment_for(distribution), distribution, *names)


def test_mismatched_nodes_and_unsupported_geometry_are_rejected():
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 1])
    experiment = experiment_for(distribution)
    reversed_distribution = ParticleSizeDistribution([300, 100] * ureg.nanometer, [1, 1])
    with pytest.raises(ValueError, match='same order'):
        average(experiment, reversed_distribution, 'Csca')
    experiment.scatterer_set = InfiniteCylinderSet(diameter=distribution.diameters, material=[1.5], medium=[1.])
    with pytest.raises(TypeError, match='homogeneous spheres'):
        average(experiment, distribution, 'Csca')


def test_lognormal_optical_quadrature_converges():
    values = []
    for sampling in (16, 32, 64):
        distribution = ParticleSizeDistribution.lognormal(100 * ureg.nanometer, 1.2, sampling)
        values.append(average(experiment_for(distribution), distribution, 'Csca').as_numpy())
    np.testing.assert_allclose(values[0], values[2], rtol=1e-7)
    np.testing.assert_allclose(values[1], values[2], rtol=1e-10)


def test_sequential_sets_are_rejected_before_averaging():
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 1])
    experiment = experiment_for(distribution)
    experiment.scatterer_set = SphereSet.build_sequential(
        diameter=distribution.diameters, material=1.5, medium=1., target_size=2,
    )
    with pytest.raises(ValueError, match='not sequential'):
        average(experiment, distribution, 'Csca')


def test_zero_scattering_has_no_defined_ensemble_asymmetry(monkeypatch):
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 1])
    experiment = experiment_for(distribution)
    # Use exactly zero scattering to exercise the mathematical degeneracy without
    # depending on floating-point residuals in the index-matched Mie recurrence.
    from PyMieSim import LabeledArray
    original = experiment._build_labeled_array

    def zero_scattering(measures, drop_unique_level):
        result = original(measures, drop_unique_level)
        if measures == ['Csca']:
            return LabeledArray(np.zeros(result.shape), list(result.dims), dict(result.coords), dict(result.attrs), result.name)
        return result

    monkeypatch.setattr(experiment, '_build_labeled_array', zero_scattering)
    with pytest.raises(ValueError, match='undefined'):
        average(experiment, distribution, 'g')


def test_detector_power_averaging_preserves_detector_axes():
    from PyMieSim import PhotodiodeSet
    distribution = ParticleSizeDistribution([100, 300] * ureg.nanometer, [1, 3])
    experiment = experiment_for(distribution)
    detector = PhotodiodeSet(
        numerical_aperture=[0.1, 0.2], sampling=[100],
        phi_offset=[0] * ureg.degree, gamma_offset=[0] * ureg.degree,
        medium=[1.0],
    )
    experiment = Experiment(source_set=experiment.source_set, scatterer_set=experiment.scatterer_set, detector_set=detector)
    raw = experiment.get('coupling')
    result = average(experiment, distribution, 'coupling', as_result=True)
    expected = np.tensordot(raw.as_numpy(), distribution.number_fractions,
                            axes=(raw.dims.index('scatterer:diameter'), 0))
    np.testing.assert_allclose(result.magnitude, expected)
    assert result.units == ureg.watt
    assert result.dims == tuple(dim for dim in raw.dims if dim != 'scatterer:diameter')


@pytest.mark.parametrize('order', [1, 2, 3, 6])
def test_uniform_moments(order):
    distribution = ParticleSizeDistribution.uniform(100 * ureg.nanometer, 0.3 * ureg.micrometer, sampling=8)
    actual = np.sum(distribution.number_fractions * distribution.diameters.to('nanometer').magnitude ** order)
    expected = (300 ** (order + 1) - 100 ** (order + 1)) / ((order + 1) * 200)
    assert actual == pytest.approx(expected, rel=1e-12)


def test_truncated_normal_matches_conditional_moments():
    from math import erf, exp, pi, sqrt
    mean, width, lower, upper = 200., 50., 150., 300.
    distribution = ParticleSizeDistribution.truncated_normal(
        mean * ureg.nanometer, width * ureg.nanometer,
        minimum_diameter=lower * ureg.nanometer, maximum_diameter=upper * ureg.nanometer,
    )
    alpha, beta = (lower - mean) / width, (upper - mean) / width
    phi_alpha = exp(-alpha ** 2 / 2) / sqrt(2 * pi)
    phi_beta = exp(-beta ** 2 / 2) / sqrt(2 * pi)
    probability = (erf(beta / sqrt(2)) - erf(alpha / sqrt(2))) / 2
    mean_shift = (phi_alpha - phi_beta) / probability
    expected_mean = mean + width * mean_shift
    expected_variance = width ** 2 * (1 + (alpha * phi_alpha - beta * phi_beta) / probability - mean_shift ** 2)
    sizes = distribution.diameters.to('nanometer').magnitude
    actual_mean = sizes @ distribution.number_fractions
    assert actual_mean == pytest.approx(expected_mean, rel=1e-12)
    assert ((sizes - actual_mean) ** 2) @ distribution.number_fractions == pytest.approx(expected_variance, rel=1e-12)
    assert sizes.min() > lower and sizes.max() < upper


@pytest.mark.parametrize('mode', [100., 140., 300.])
def test_triangular_moments_include_endpoint_modes(mode):
    lower, upper = 100., 300.
    distribution = ParticleSizeDistribution.triangular(
        lower * ureg.nanometer, mode * ureg.nanometer, upper * ureg.nanometer, sampling=8,
    )
    sizes = distribution.diameters.to('nanometer').magnitude
    actual_mean = sizes @ distribution.number_fractions
    expected_variance = (lower ** 2 + upper ** 2 + mode ** 2 - lower * upper - lower * mode - upper * mode) / 18
    assert actual_mean == pytest.approx((lower + mode + upper) / 3)
    assert (sizes - actual_mean) ** 2 @ distribution.number_fractions == pytest.approx(expected_variance)


def test_mixture_uses_component_number_fractions_not_node_counts():
    small = ParticleSizeDistribution.monodisperse(100 * ureg.nanometer)
    large = ParticleSizeDistribution.uniform(200 * ureg.nanometer, 400 * ureg.nanometer, sampling=16)
    mixture = ParticleSizeDistribution.mixture([small, large], [3, 1])
    assert mixture.number_fractions[0] == pytest.approx(0.75)
    assert mixture.number_fractions[1:].sum() == pytest.approx(0.25)
    assert mixture.diameters.to('nanometer').magnitude @ mixture.number_fractions == pytest.approx(150)
    result = average(experiment_for(mixture), mixture, 'Csca').as_numpy()
    expected = (0.75 * average(experiment_for(small), small, 'Csca').as_numpy()
                + 0.25 * average(experiment_for(large), large, 'Csca').as_numpy())
    np.testing.assert_allclose(result, expected, rtol=1e-12)


@pytest.mark.parametrize('family', ['uniform', 'truncated_normal', 'triangular'])
def test_new_continuous_families_converge_in_optical_averages(family):
    def distribution(sampling):
        bounds = {'minimum_diameter': 100 * ureg.nanometer, 'maximum_diameter': 300 * ureg.nanometer}
        if family == 'truncated_normal':
            bounds.update(mean_diameter=200 * ureg.nanometer, standard_deviation=40 * ureg.nanometer)
        elif family == 'triangular':
            bounds.update(mode_diameter=180 * ureg.nanometer)
        return getattr(ParticleSizeDistribution, family)(**bounds, sampling=sampling)
    coarse, fine = distribution(16), distribution(32)
    np.testing.assert_allclose(average(experiment_for(coarse), coarse, 'Csca').as_numpy(),
                               average(experiment_for(fine), fine, 'Csca').as_numpy(), rtol=1e-8)


@pytest.mark.parametrize('family', ['uniform', 'truncated_normal', 'triangular'])
@pytest.mark.parametrize('invalid', [{'minimum_diameter': 0 * ureg.nanometer},
                                    {'maximum_diameter': 50 * ureg.nanometer},
                                    {'sampling': 0}, {'sampling': True}, {'sampling': 2.5}])
def test_invalid_continuous_family_parameters(family, invalid):
    kwargs = {'minimum_diameter': 100 * ureg.nanometer, 'maximum_diameter': 300 * ureg.nanometer}
    if family == 'truncated_normal':
        kwargs.update(mean_diameter=200 * ureg.nanometer, standard_deviation=40 * ureg.nanometer)
    elif family == 'triangular':
        kwargs.update(mode_diameter=180 * ureg.nanometer)
    with pytest.raises(ValueError):
        getattr(ParticleSizeDistribution, family)(**(kwargs | invalid))


def test_normal_requires_positive_width_and_explicit_bounds():
    with pytest.raises(TypeError):
        ParticleSizeDistribution.truncated_normal(200 * ureg.nanometer, 40 * ureg.nanometer)
    with pytest.raises(ValueError, match='standard_deviation'):
        ParticleSizeDistribution.truncated_normal(200 * ureg.nanometer, 0 * ureg.nanometer,
                                                  minimum_diameter=100 * ureg.nanometer, maximum_diameter=300 * ureg.nanometer)
    with pytest.raises(ValueError, match='mode_diameter'):
        ParticleSizeDistribution.triangular(100 * ureg.nanometer, 400 * ureg.nanometer, 300 * ureg.nanometer)


@pytest.mark.parametrize('components, weights', [([], []), ([None], [1]),
    ([ParticleSizeDistribution.monodisperse(100 * ureg.nanometer)], [0]),
    ([ParticleSizeDistribution.monodisperse(100 * ureg.nanometer)], [1, 1])])
def test_invalid_mixture_components_and_weights(components, weights):
    with pytest.raises(ValueError):
        ParticleSizeDistribution.mixture(components, weights)
