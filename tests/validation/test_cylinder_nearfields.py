import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from PyMieSim.polarization import PolarizationState
from PyMieSim.single import Setup
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import ureg


def _setup(material=1.5, angle=0.0):
    source = PlaneWave(
        wavelength=600 * ureg.nanometer,
        polarization=PolarizationState(angle=angle * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    scatterer = InfiniteCylinder(
        diameter=300 * ureg.nanometer,
        material=material,
        medium=1.0,
    )
    return Setup(scatterer=scatterer, source=source)


def _coordinates(radius):
    x = np.array([radius]) * ureg.nanometer
    y = np.array([0.0]) * ureg.nanometer
    z = np.array([0.0]) * ureg.nanometer
    return x, y, z


def test_cylinder_nearfields_reduce_to_incident_field_without_contrast():
    setup = _setup(material=1.0, angle=0.0)
    x, y, z = _coordinates(100.0)

    total_ex = setup.get_total_nearfields(x, y, z, "Ex").magnitude
    total_ez = setup.get_total_nearfields(x, y, z, "Ez").magnitude
    scattered_ex = setup.get_scattered_nearfields(x, y, z, "Ex").magnitude

    assert np.allclose(total_ex, 1.0, atol=1e-12)
    assert np.allclose(total_ez, 0.0, atol=1e-12)
    assert np.allclose(scattered_ex, 0.0, atol=1e-12)


def test_cylinder_tangential_fields_are_continuous_at_surface():
    radius = 150.0

    setup_x = _setup(angle=0.0)
    x_inside, y_inside, z_inside = _coordinates(radius - 1e-3)
    x_outside, y_outside, z_outside = _coordinates(radius + 1e-3)

    ez_inside = setup_x.get_total_nearfields(x_inside, y_inside, z_inside, "Ez").magnitude
    ez_outside = setup_x.get_total_nearfields(x_outside, y_outside, z_outside, "Ez").magnitude
    assert np.allclose(ez_inside, ez_outside, rtol=1e-4, atol=1e-6)

    setup_y = _setup(angle=90.0)
    ey_inside = setup_y.get_total_nearfields(x_inside, y_inside, z_inside, "Ey").magnitude
    ey_outside = setup_y.get_total_nearfields(x_outside, y_outside, z_outside, "Ey").magnitude
    assert np.allclose(ey_inside, ey_outside, rtol=1e-4, atol=1e-6)


def test_cylinder_exterior_total_field_decomposes():
    setup = _setup(angle=0.0)
    x = np.array([0.4]) * ureg.micrometer
    y = np.array([0.0]) * ureg.micrometer
    z = np.array([0.1]) * ureg.micrometer

    total = setup.get_total_nearfields(x, y, z, "Ex").magnitude
    scattered = setup.get_scattered_nearfields(x, y, z, "Ex").magnitude
    incident = setup.get_incident_nearfields(x, y, z, "Ex").magnitude

    assert np.allclose(total, incident + scattered, rtol=1e-10, atol=1e-12)


def test_cylinder_outline_matches_plane_projection():
    nearfields = _setup().get_representation("nearfields")

    nearfields._setup_coordinates_plane(
        sampling=4,
        plane_normal=(1.0, 1.0, 0.0),
    )
    figure, axis = plt.subplots()
    nearfields._add_scatterer_outline_on_plane(axis)
    assert len(axis.patches) == 1
    assert isinstance(axis.patches[0], Ellipse)
    plt.close(figure)

    nearfields._setup_coordinates_plane(
        sampling=4,
        plane_normal=(1.0, 0.0, 0.0),
    )
    figure, axis = plt.subplots()
    nearfields._add_scatterer_outline_on_plane(axis)
    assert len(axis.lines) == 2
    plt.close(figure)
