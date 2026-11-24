"""Check the polarization, phase, and absolute power of sphere far fields."""

if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__]))

import numpy as np
import pytest
from PyMieSim import Experiment
from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.polarization import PolarizationState
from PyMieSim.single import Simulation
from PyMieSim.single.detector import CoherentMode, IntegratingSphere, Photodiode
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.setup import Setup
from PyMieSim.single.source import Gaussian, PlaneWave
from PyMieSim.units import ureg


SPEED_OF_LIGHT = 299_792_458.0  # m/s
VACUUM_PERMITTIVITY = 8.854_187_818_8e-12  # F/m


@pytest.mark.parametrize("structured", [False, True])
def test_farfields_use_parallel_and_perpendicular_mie_amplitudes(structured):
    """An x-polarized wave selects different complex amplitudes in two azimuths."""
    wavelength = 600 * ureg.nanometer
    amplitude = 2.3 * ureg.volt / ureg.meter
    distance = 1 * ureg.millimeter
    medium_index = 1.33
    sphere = Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=medium_index)
    source = PlaneWave(wavelength=wavelength, polarization=PolarizationState(angle=0 * ureg.degree), amplitude=amplitude)
    simulation = Simulation(scatterer=sphere, source=source)
    perpendicular_amplitude, parallel_amplitude = (
        values.magnitude[0] for values in simulation.get_s1s2(angles=[60] * ureg.degree)
    )

    if structured:
        # On the 13-point mesh, row 8 is elevation 30°, and columns 6 and 9 are azimuths 0° and 90°.
        elevation_field, azimuth_field, _ = simulation.get_farfields(sampling=13, distance=distance)
        elevation_field = elevation_field[8, [6, 9]]
        azimuth_field = azimuth_field[8, [6, 9]]
    else:
        elevation_field, azimuth_field = simulation.get_farfields(
            phi=[30, 30] * ureg.degree, theta=[0, 90] * ureg.degree, distance=distance
        )

    elevation = elevation_field.to("volt/meter").magnitude
    azimuth = azimuth_field.to("volt/meter").magnitude
    field_scale = amplitude.magnitude / (2 * np.pi * medium_index * (distance / wavelength).to("").magnitude)
    np.testing.assert_allclose(
        [abs(elevation[0]), abs(azimuth[1])],
        field_scale * np.abs([parallel_amplitude, perpendicular_amplitude]),
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose([azimuth[0], elevation[1]], 0, atol=1e-12 * field_scale)
    np.testing.assert_allclose(elevation[0] / azimuth[1], parallel_amplitude / perpendicular_amplitude, rtol=1e-12)


@pytest.mark.parametrize("medium_index", [1.0, 1.33])
def test_farfield_advances_by_quarter_cycle_in_host(medium_index):
    """A quarter of the host wavelength rotates the radius-scaled field by +i."""
    wavelength = 600 * ureg.nanometer
    simulation = Simulation(
        scatterer=Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=medium_index),
        source=PlaneWave(
            wavelength=wavelength,
            polarization=PolarizationState(angle=0 * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        ),
    )
    first_radius = 100 * ureg.micrometer
    second_radius = first_radius + wavelength / (4 * medium_index)
    scaled_fields = []
    for radius in (first_radius, second_radius):
        elevation_field, _ = simulation.get_farfields(phi=[30] * ureg.degree, theta=[0] * ureg.degree, distance=radius)
        scaled_fields.append((radius * elevation_field[0]).to("volt").magnitude)

    np.testing.assert_allclose(scaled_fields[1] / scaled_fields[0], 1j, rtol=0, atol=1e-11)


def test_sphere_farfield_matches_complex_nearfield_at_large_radius():
    """The scattered-field asymptote agrees in both transverse components."""
    setup = Setup(
        scatterer=Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=1.33),
        source=PlaneWave(
            wavelength=600 * ureg.nanometer,
            polarization=PolarizationState(angle=35 * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        ),
    )
    elevation = np.pi / 6
    radius = 1 * ureg.millimeter
    coordinates = dict(
        x=np.array([np.cos(elevation)]) * radius,
        y=np.array([0.0]) * radius,
        z=np.array([np.sin(elevation)]) * radius,
    )
    nearfield_x, nearfield_y, nearfield_z = (
        setup.get_scattered_nearfields(**coordinates, field_type=component).to("volt/meter").magnitude[0]
        for component in ("Ex", "Ey", "Ez")
    )
    # At azimuth zero, the elevation unit vector is (-sin(elevation), 0, cos(elevation))
    # and the azimuth unit vector is (0, 1, 0).
    nearfield_transverse = [
        -np.sin(elevation) * nearfield_x + np.cos(elevation) * nearfield_z,
        nearfield_y,
    ]
    farfield_transverse = [
        component.to("volt/meter").magnitude[0]
        for component in setup.get_farfields(phi=[elevation] * ureg.radian, theta=[0] * ureg.radian, distance=radius)
    ]
    # The first neglected radial term is of relative order 1/(k_medium * radius).
    np.testing.assert_allclose(farfield_transverse, nearfield_transverse, rtol=5e-4, atol=0)


@pytest.mark.parametrize("medium_index", [1.0, 1.33])
def test_integrating_sphere_power_equals_irradiance_times_cross_section(medium_index):
    """Full-sphere power follows the scattering cross section and peak incident field."""
    amplitude = 2.3 * ureg.volt / ureg.meter
    sphere = Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=medium_index)
    source = PlaneWave(
        wavelength=600 * ureg.nanometer, polarization=PolarizationState(angle=0 * ureg.degree), amplitude=amplitude
    )
    detector = IntegratingSphere(sampling=20_000)
    simulation = Simulation(scatterer=sphere, source=source, detector=detector)

    cross_section = simulation.run("Csca").quantity.to("meter**2").magnitude
    irradiance = 0.5 * medium_index * VACUUM_PERMITTIVITY * SPEED_OF_LIGHT * amplitude.magnitude**2
    expected_power = irradiance * cross_section
    # The Fibonacci mesh approximates the full-sphere integral at finite sampling.
    for collected_power in (
        simulation.run("coupling").quantity,
        detector.get_energy_flow(scatterer=sphere, source=source),
    ):
        np.testing.assert_allclose(collected_power.to("watt").magnitude, expected_power, rtol=5e-4, atol=0)


@pytest.mark.parametrize("polarization_angle", [0, 90])
def test_experiment_photodiode_power_matches_differential_cross_section(polarization_angle):
    """A narrow side detector resolves polarized power in the experiment API."""
    wavelength = 600 * ureg.nanometer
    diameter = 400 * ureg.nanometer
    angle = polarization_angle * ureg.degree
    optical_power = 1 * ureg.watt
    source_aperture = 0.3
    detector_aperture = 0.005
    simulation = Simulation(
        scatterer=Sphere(diameter=diameter, material=1.6, medium=1.0),
        source=Gaussian(
            wavelength=wavelength,
            polarization=PolarizationState(angle=angle),
            optical_power=optical_power,
            numerical_aperture=source_aperture,
        ),
    )
    perpendicular_amplitude, parallel_amplitude = (
        values.magnitude[0] for values in simulation.get_s1s2(angles=[90] * ureg.degree)
    )
    # The detector is in the x-z plane, so x polarization selects S2 and y selects S1.
    scattering_amplitude = parallel_amplitude if polarization_angle == 0 else perpendicular_amplitude

    source_set = GaussianSet(
        wavelength=wavelength,
        polarization=PolarizationSet(angles=angle),
        optical_power=optical_power,
        numerical_aperture=source_aperture,
    )
    experiment = Experiment(
        source_set=source_set,
        scatterer_set=SphereSet(diameter=diameter, material=1.6, medium=1.0),
        detector_set=PhotodiodeSet(
            numerical_aperture=detector_aperture, phi_offset=90 * ureg.degree, gamma_offset=0 * ureg.degree, sampling=2000
        ),
    )
    collected_power = experiment.run("coupling").quantity.to("watt").magnitude.item()
    vacuum_wavenumber = 2 * np.pi / wavelength.to("meter").magnitude
    beam_waist = 0.61 * wavelength.to("meter").magnitude / source_aperture
    peak_irradiance = 2 * optical_power.magnitude / (np.pi * beam_waist**2)
    cone_solid_angle = 2 * np.pi * (1 - np.sqrt(1 - detector_aperture**2))
    expected_power = peak_irradiance * abs(scattering_amplitude / vacuum_wavenumber)**2 * cone_solid_angle
    # Use the cone center as the differential cross-section reference.
    np.testing.assert_allclose(collected_power, expected_power, rtol=1e-3, atol=0)


def test_side_photodiode_projects_parallel_field_onto_its_rotated_axis():
    """The analyzer along the rotated x axis selects in-plane polarization."""
    source = PlaneWave(
        wavelength=600 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    powers = []
    for filter_angle in (0, 90):
        detector = Photodiode(
            numerical_aperture=0.005,
            phi_offset=90 * ureg.degree,
            gamma_offset=0 * ureg.degree,
            polarization_filter=filter_angle * ureg.degree,
            sampling=2000,
        )
        simulation = Simulation(
            scatterer=Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=1.0),
            source=source,
            detector=detector,
        )
        powers.append(simulation.run("coupling").quantity.to("watt").magnitude)

    parallel_power, perpendicular_power = powers
    assert parallel_power > 0
    # A finite cone admits a small cross-polarized component off its central ray.
    assert perpendicular_power / parallel_power < 1e-4


@pytest.mark.parametrize("circular", [False, True])
def test_forward_photodiode_analyzer_transmits_matching_polarization(circular):
    """Aligned and orthogonal analyzers distinguish relative phase in the field."""
    if circular:
        weight = 1 / np.sqrt(2)
        incident = PolarizationState(jones_vector=[weight, 1j * weight])
        crossed = PolarizationState(jones_vector=[weight, -1j * weight])
    else:
        incident = PolarizationState(angle=45 * ureg.degree)
        crossed = PolarizationState(angle=135 * ureg.degree)

    source = PlaneWave(wavelength=600 * ureg.nanometer, polarization=incident, amplitude=1 * ureg.volt / ureg.meter)
    powers = []
    for analyzer in (None, incident, crossed):
        filter_options = {} if analyzer is None else {"polarization_filter": analyzer}
        detector = Photodiode(
            numerical_aperture=0.001,
            phi_offset=0 * ureg.degree,
            gamma_offset=0 * ureg.degree,
            sampling=2000,
            **filter_options,
        )
        simulation = Simulation(
            scatterer=Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=1.0),
            source=source,
            detector=detector,
        )
        powers.append(simulation.run("coupling").quantity.to("watt").magnitude)

    unfiltered_power, aligned_power, crossed_power = powers
    assert unfiltered_power > 0
    # A finite cone admits a small cross-polarized contribution away from its axis.
    np.testing.assert_allclose(aligned_power / unfiltered_power, 1, rtol=0, atol=1e-6)
    assert crossed_power / unfiltered_power < 1e-6


def test_experiment_photodiode_sweep_resolves_crossed_analyzers():
    """An analyzer sweep follows Malus's law for a narrow forward cone."""
    experiment = Experiment(
        scatterer_set=SphereSet(diameter=400 * ureg.nanometer, material=1.6, medium=1.0),
        source_set=PlaneWaveSet(
            wavelength=600 * ureg.nanometer,
            polarization=PolarizationSet(angles=[45] * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        ),
        detector_set=PhotodiodeSet(
            numerical_aperture=0.001,
            phi_offset=0 * ureg.degree,
            gamma_offset=0 * ureg.degree,
            polarization_filter=[45, 90, 135] * ureg.degree,
            sampling=2000,
        ),
    )
    powers = experiment.run("coupling").as_numpy()
    assert powers[0] > 0
    np.testing.assert_allclose(powers / powers[0], [1, 0.5, 0], rtol=0, atol=1e-6)


@pytest.mark.parametrize("mean_coupling", [False, True])
def test_coherent_mode_power_uses_host_irradiance(mean_coupling):
    """The mode measures peak-phasor power in a non-vacuum surrounding medium."""
    wavelength = 600 * ureg.nanometer
    amplitude = 1 * ureg.volt / ureg.meter
    medium_index = 1.33
    cone_half_angle = 0.005
    sampling = 2000
    detector = CoherentMode(
        mode_number="LP01",
        numerical_aperture=medium_index * np.sin(cone_half_angle),
        phi_offset=60 * ureg.degree,
        gamma_offset=0 * ureg.degree,
        rotation=0 * ureg.degree,
        medium=medium_index,
        mean_coupling=mean_coupling,
        sampling=sampling,
    )
    simulation = Simulation(
        scatterer=Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=medium_index),
        source=PlaneWave(
            wavelength=wavelength, polarization=PolarizationState(angle=0 * ureg.degree), amplitude=amplitude
        ),
        detector=detector,
    )
    _, parallel_amplitude = simulation.get_s1s2(angles=[60] * ureg.degree)
    medium_wavenumber = 2 * np.pi * medium_index / wavelength.to("meter").magnitude
    irradiance = 0.5 * medium_index * VACUUM_PERMITTIVITY * SPEED_OF_LIGHT * amplitude.magnitude**2
    power_per_solid_angle = irradiance * abs(parallel_amplitude.magnitude[0] / medium_wavenumber)**2
    solid_angle_per_sample = 2 * np.pi * (1 - np.cos(cone_half_angle)) / sampling
    mode = detector.scalar_field
    if mean_coupling:
        expected_power = power_per_solid_angle * np.sum(np.abs(mode)**2) / sampling
    else:
        expected_power = power_per_solid_angle * abs(np.sum(mode))**2 * solid_angle_per_sample
    collected_power = simulation.run("coupling").quantity.to("watt").magnitude
    np.testing.assert_allclose(collected_power, expected_power, rtol=1e-3, atol=0)
