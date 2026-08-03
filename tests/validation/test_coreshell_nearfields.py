import numpy as np
import matplotlib.pyplot as plt

from PyMieSim.polarization import PolarizationState
from PyMieSim.single import Setup
from PyMieSim.single.representations import NearFields
from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import ureg


def _setup():
    source = PlaneWave(
        wavelength=600 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    scatterer = CoreShell(
        core_diameter=200 * ureg.nanometer,
        shell_thickness=100 * ureg.nanometer,
        core_material=1.5,
        shell_material=2.0,
        medium=1.0,
    )
    return Setup(scatterer=scatterer, source=source)


def test_coreshell_nearfields_are_finite_in_all_radial_regions():
    setup = _setup()
    coordinates = np.array([0.05, 0.15, 0.25, 0.4]) * ureg.micrometer
    zeros = np.zeros(4) * ureg.micrometer

    total = setup.get_total_nearfields(coordinates, zeros, zeros, "Ex").magnitude
    scattered = setup.get_scattered_nearfields(coordinates, zeros, zeros, "Ex").magnitude

    assert np.all(np.isfinite(total))
    assert np.all(np.isfinite(scattered))
    assert np.any(np.abs(total[:3]) > 0)
    assert np.allclose(scattered[:2], 0.0)


def test_coreshell_exterior_total_field_decomposes_into_incident_and_scattered():
    setup = _setup()
    x = np.array([0.4, 0.5]) * ureg.micrometer
    y = np.array([0.1, 0.2]) * ureg.micrometer
    z = np.array([0.2, 0.1]) * ureg.micrometer

    incident = setup.get_incident_nearfields(x, y, z, "Ex").magnitude
    total = setup.get_total_nearfields(x, y, z, "Ex").magnitude
    scattered = setup.get_scattered_nearfields(x, y, z, "Ex").magnitude

    assert np.allclose(total, incident + scattered, rtol=1e-10, atol=1e-12)


def test_coreshell_nearfield_representation_infers_extent_from_total_diameter():
    nearfields = NearFields(_setup())

    u_range, v_range = nearfields._default_uv_range(extent_scale=2.5)

    assert u_range[0].to(ureg.nanometer).magnitude == -500
    assert u_range[1].to(ureg.nanometer).magnitude == 500
    assert v_range == u_range


def test_coreshell_nearfield_plot_draws_core_and_outer_shell_circles():
    nearfields = NearFields(_setup())
    nearfields._setup_coordinates_plane(sampling=4)
    figure, axis = plt.subplots()

    nearfields._add_scatterer_outline_on_plane(axis)

    assert len(axis.patches) == 2
    assert {patch.get_linestyle() for patch in axis.patches} == {"--", "-"}
    plt.close(figure)
