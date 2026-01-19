import pytest
import numpy as np
from PyMieSim.units import ureg, Angle
from PyMieSim.single.mesh import FibonacciMesh
from PyMieSim.units import ureg



mesh = FibonacciMesh(
    max_angle=np.pi / 4 * ureg.radian,  # Maximum angle (45 degrees in radians)
    sampling=100 * ureg.AU,  # Number of sample points
    phi_offset=30 * ureg.degree,  # Offset in the azimuthal angle (in degrees)
    gamma_offset=45 * ureg.degree,  # Offset in the polar angle (in degrees)
    rotation_angle=0 * ureg.degree,  # Rotation angle (in degrees)
)