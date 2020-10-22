
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

class S1S2Plot(object):

    def __init__(self, S1, S2, X, Y, Z, Mesh):

        self.S2 = S2
        self.S1 = S1
        self.Mesh = Mesh

        self.fig = plt.figure(figsize=(15, 5))

        self.ax0 = self.fig.add_subplot(131, projection='polar')
        self.ax1 = self.fig.add_subplot(132, projection='polar')
        self.ax2 = self.fig.add_subplot(133, projection='3d')

        self.ax0.set_title(r'$S_1$ [logscale]')
        self.ax1.set_title(r'$S_2$ [logscale]')
        self.ax2.set_title(r'3D Phase function')


        data = np.abs(self.S1)
        self.ax0.plot(self.Mesh.PhiVec.Radian,
                 data,
                 'k')

        self.ax0.fill_between(self.Mesh.PhiVec.Radian,
                         0,
                         data,
                         color='C0',
                         alpha=0.4)



        data = np.abs(self.S2)
        self.ax1.plot(self.Mesh.PhiVec.Radian,
                 data,
                 'k')

        self.ax1.fill_between(self.Mesh.PhiVec.Radian,
                         0,
                         data,
                         color='C1',
                         alpha=0.4)

        self.ax2.plot_surface(Y,
                         X,
                         Z,
                         rstride=3,
                         cstride=3,
                         linewidth=0.5,
                         cmap=cm.bone,
                         antialiased=False,
                         alpha=1)

        self.ax2.set_xlabel('X direction')
        self.ax2.set_ylabel('Y direction')
        self.ax2.set_zlabel('Z direction')
        norm = np.sqrt(np.max(X**2+Y**2))

        self.ax2.set_xlim([-norm, norm])
        self.ax2.set_ylim([-norm, norm])



        plt.show()
