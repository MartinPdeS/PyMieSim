import matplotlib.pyplot as plt
from PyMieSim.single.mesh import FibonacciMesh
from TypedUnit import ureg


mesh = FibonacciMesh(
    sampling=300,
    min_angle=30 * ureg.degree,
    max_angle=50 * ureg.degree,
    phi_offset=0 * ureg.degree,
    gamma_offset=0 * ureg.degree,
)

print(mesh.spherical.phi.shape)


figure = plt.figure()

ax = figure.add_subplot(111, projection="3d")


ax.scatter(mesh.cartesian.x, mesh.cartesian.y, mesh.cartesian.z)

max_abs = 1
ax.set(xlim=(-max_abs, max_abs), ylim=(-max_abs, max_abs), zlim=(-max_abs, max_abs))
plt.show()
