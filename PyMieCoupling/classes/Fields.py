
from numpy import abs, arctan, angle, sqrt, real, imag, conjugate, exp, array, max
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patches as mpatches

class Field(object):

    def __init__(self, Perpendicular, Parallel, Meshes):


        self.__dict__ = Meshes.__dict__.copy()

        self.Perpendicular = Perpendicular

        self.Parallel = Parallel

        self.delta = angle(self.Parallel)-angle(self.Perpendicular)

        self.Meshes = Meshes

        self.SPF = abs(Parallel)**2 + abs(Perpendicular)**2

        ## Source -- https://www.physlab.org/wp-content/uploads/2016/07/Ch6-BYUOpticsBook_2013.pdf

        self.Total = sqrt(abs(Parallel)**2 + abs(Perpendicular)**2)# * exp(complex(0,1)*self.delta)

        self.Angle = arctan(abs(Parallel)/abs(Perpendicular))

        self.Jones, self.Stokes = None, None

        self.ComputeJones()

        self.ComputeStokes()



    def ComputeStokes(self):

        I = abs(self.Parallel)**2 + abs(self.Perpendicular)**2

        self.Stokes = array( [ I,
                               (abs(self.Parallel)**2 - abs(self.Perpendicular)**2)/I,
                                2*real( self.Parallel * conjugate(self.Perpendicular))/I,
                                -2*imag( self.Parallel * conjugate(self.Perpendicular))/I
                                ])



    def ComputeJones(self):

        A = abs(self.Parallel) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        B = abs(self.Perpendicular) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        self.Jones = array([A, B*exp(complex(0,1)*self.delta)])



    def GenStokesFig(self):
        fig = plt.figure(figsize=(15, 6))

        ax0 = fig.add_subplot(121)
        ax1 = fig.add_subplot(122)

        ax0.set_title(r'$S_0$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$] ')
        ax0.set_xlabel(r'Angle $\phi$')
        ax0.set_ylabel(r'Angle $\theta$')

        ax1.set_title(r'$S_3$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$]')
        ax1.set_ylabel(r'Angle $\theta$')

        return fig, ax0, ax1


    def PlotStokes(self, RectangleTheta=[0,360], RectanglePhi=[0,360]):

        Origin = (RectangleTheta[1], RectanglePhi[1])
        SolidAngle = (RectangleTheta[0] - RectangleTheta[1], RectanglePhi[0] - RectanglePhi[1])

        RectangleTheta = array(RectangleTheta) % 360


        fig, ax0, ax1 = self.GenStokesFig()

        im0 = ax0.pcolormesh(
                             self.Phi.Vector.Degree,
                             self.Theta.Vector.Degree,
                             self.Stokes[0].T,
                             shading='auto'
                             )

        cbar = plt.colorbar(im0, ax=ax0, pad=0.15, orientation='horizontal')
        cbar.ax.tick_params(labelsize='small')
        cbar.ax.locator_params(nbins=3)

        n=8
        im1 = ax0.quiver(
                         self.Phi.Mesh.Degree[::n, ::n],
                         self.Theta.Mesh.Degree[::n, ::n],
                         self.Stokes[2][::n, ::n],
                         self.Stokes[1][::n, ::n],
                         units='width',
                         width=0.0005,
                         headwidth=30,
                         headlength=20,
                         headaxislength=20,
                         )


        im2 = ax1.pcolormesh(
                             self.Phi.Vector.Degree,
                             self.Theta.Vector.Degree,
                             self.Stokes[3].T,
                             shading='auto'
                             )

        cbar = plt.colorbar(im2, ax=ax1, pad=0.15, orientation='horizontal')
        cbar.ax.tick_params(labelsize='small')
        cbar.ax.locator_params(nbins=3)

        ax1.quiver(
                   self.Phi.Mesh.Degree[::n, ::n],
                   self.Theta.Mesh.Degree[::n, ::n],
                   self.Stokes[2][::n, ::n],
                   self.Stokes[1][::n, ::n],
                   units='width',
                   width=0.0005,
                   headwidth=30,
                   headlength=20,
                   headaxislength=20,
                         )


        rect0 = patches.Rectangle(Origin,*SolidAngle,linewidth=1,edgecolor='r',facecolor='none')

        rect1 = patches.Rectangle(Origin,*SolidAngle,linewidth=1,edgecolor='r',facecolor='none')

        ax0.add_patch(rect0)

        ax1.add_patch(rect1)


        plt.show()














# -
