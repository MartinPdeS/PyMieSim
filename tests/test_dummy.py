from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian


def test_simple():
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        index=1.5,
        n_medium=1.0
    )

    # print(scatterer.Qsca)

# -
