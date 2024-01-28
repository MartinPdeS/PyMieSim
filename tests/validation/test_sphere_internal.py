import numpy

from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.mesh import FibonacciMesh


def test_Qsca_cross_section():
    source = Gaussian(
        wavelength=1e-6,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )
    sphere = Sphere(diameter=300e-9, index=1.4, source=source)
    val0 = sphere.get_cross_section()
    val1 = sphere.Qsca * sphere.area

    if not numpy.isclose(val0, val1, atol=0, rtol=1e-5):
        raise ValueError('Mismatch with testing values')


def test_energy_flow_coupling():
    source = Gaussian(
        wavelength=1e-6,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )
    sphere = Sphere(
        diameter=300e-9,
        index=1.4,
        source=source
    )

    detector = Photodiode(
        sampling=500,
        NA=2.0,
        gamma_offset=0,
        phi_offset=0
    )

    val0 = sphere.get_energy_flow(detector.cpp_binding.mesh)
    val1 = detector.coupling(sphere)

    if not numpy.isclose(val0, val1, atol=0, rtol=1e-5):
        raise ValueError('Mismatch with testing values')


def test_compare_sphere_coreshell_0():
    source = Gaussian(
        wavelength=1e-6,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    sphere = Sphere(
        diameter=1e-6,
        index=1.5,
        source=source
    )

    coreshell = CoreShell(
        core_diameter=1e-6,
        shell_width=0,
        core_index=1.5,
        shell_index=1.8,
        source=source
    )

    for parameter in ['Qsca', 'Qext', 'Qabs']:
        value_sphere = getattr(sphere, parameter)
        value_coreshell = getattr(coreshell, parameter)

        if not numpy.isclose(value_sphere, value_coreshell, atol=1e-12, rtol=1e-5):
            raise ValueError(f'Mismatch with juxtaposing CoreShell with zero shell and Sphere for parameter: {parameter}')

# -
