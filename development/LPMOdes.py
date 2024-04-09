from PyMieSim.single.detector import HGMode

test = HGMode(
    mode_number='HG00',
    NA=0.2,
    gamma_offset=0,
    phi_offset=0
)

# from PyMieSim.binary.DetectorInterface import BindedDetector


# from MPSPlots.render3D import SceneList as SceneList3D
# from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh
# import numpy

# from scipy.linalg import norm


# from PyMieSim.tools.modes.laguerre_gauss import interpolate_from_fibonacci_mesh


# fibonacci_mesh = CPPFibonacciMesh(
#     sampling=1500,
#     max_angle=0.3,
#     phi_offset=0,
#     gamma_offset=0,
#     rotation_angle=0
# )


# coordinate = numpy.row_stack((
#     fibonacci_mesh.x,
#     fibonacci_mesh.y,
#     fibonacci_mesh.z
# ))

# mode_field = interpolate_from_fibonacci_mesh(
#     fibonacci_mesh=fibonacci_mesh,
#     p=3,
#     l=0,
#     waist_radius=0.08,
# )
# print("Mode Field Amplitudes:", norm(mode_field))

# figure = SceneList3D()


# ax = figure.append_ax()

# ax.add_unstructured_mesh(
#     coordinates=coordinate,
#     scalar_coloring=mode_field,
#     symmetric_map=True,
#     symmetric_colormap=True
# )

# figure.show()


# -
