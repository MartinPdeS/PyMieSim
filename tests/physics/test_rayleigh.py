"""Small dielectric spheres approach the electric-dipole scattering limit."""

import numpy as np
from PyMieSim import Experiment
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import PlaneWaveSet
from PyMieSim.units import ureg


def test_rayleigh_scattering_efficiency_in_dielectric_host():
    """Check the leading coefficient and fourth-power diameter dependence."""
    diameters = np.array([1., 2.])
    host, particle = 1.33, 1.6
    wavelength = 600.
    source = PlaneWaveSet(
        wavelength=wavelength * ureg.nanometer, polarization=PolarizationSet(angles=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    spheres = SphereSet(diameter=diameters * ureg.nanometer, material=particle, medium=host)
    efficiency = Experiment(source_set=source, scatterer_set=spheres).run("Qsca").as_numpy()
    size = np.pi * host * diameters / wavelength
    relative_permittivity = (particle / host)**2
    expected = (8 / 3) * size**4 * ((relative_permittivity - 1) / (relative_permittivity + 2))**2
    # The first relative correction is of order x² (below 2e-4 here).
    np.testing.assert_allclose(efficiency, expected, rtol=2e-4, atol=0)
    np.testing.assert_allclose(efficiency[1] / efficiency[0], 16, rtol=2e-4, atol=0)
