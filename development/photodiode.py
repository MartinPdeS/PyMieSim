from TypedUnit import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import plot_system

detector = CoherentMode(
    mode_number="LG23",  # Specifying LP23 mode
    sampling=400 * ureg.AU,  # Number of sampling points
    NA=0.4 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset
    phi_offset=40 * ureg.degree,  # Phi offset in degrees
)

# print(detector._cpp_mesh.spherical.phi.shape)


detector._cpp_scalar_field[:200] *= 0


field = detector.get_structured_scalarfield(sampling=10)
print(field)

plot_system(detector)
print("done")
