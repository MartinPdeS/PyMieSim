
from numpy import abs, arctan, angle, sqrt, real, imag, conjugate, exp, array, max

import matplotlib.pyplot as plt

class Field(object):

    def __init__(self, Perpendicular, Parallel, Meshes):
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

        I = max(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        self.Stokes = array([abs(self.Parallel)**2 + abs(self.Perpendicular)**2/I,
                                abs(self.Parallel)**2 - abs(self.Perpendicular)**2/I,
                                2 * real( self.Parallel * conjugate(self.Perpendicular) )/I,
                                -2 * imag( self.Parallel * conjugate(self.Perpendicular) )/I
                                ])



    def ComputeJones(self):

        A = abs(self.Parallel) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        B = abs(self.Perpendicular) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        self.Jones = array([A, B*exp(complex(0,1)*self.delta)])



    def GenStokesFig(self):
        fig = plt.figure(figsize=(15, 6))

        ax0 = fig.add_subplot(121)
        ax1 = fig.add_subplot(122)

        ax0.set_title(r'$S_1$ Stokes parameter of Far-Field & Projection on [$S_2, S_3$] ')
        ax0.set_xlabel(r'Angle $\theta$')
        ax0.set_ylabel(r'Angle $\phi$')

        ax1.set_title(r'$S_4$ Stokes parameter of Far-Field & Projection on [$S_2, S_3$]')
        ax1.set_ylabel(r'Angle $\phi$')

        return fig, ax0, ax1


    def PlotStokes(self):

        fig, ax0, ax1 = self.GenStokesFig()

        im0 = ax0.pcolormesh(self.Meshes.ThetaVec.Degree,
                             self.Meshes.PhiVec.Degree,
                             self.Stokes[0],
                             shading='auto')

        cbar = plt.colorbar(im0, ax=ax0, pad=0.15, orientation='horizontal')
        cbar.ax.tick_params(labelsize='small')
        cbar.ax.locator_params(nbins=3)


        im1 = ax0.quiver(self.Meshes.ThetaMesh.Degree,
                         self.Meshes.PhiMesh.Degree,
                         self.Stokes[1],
                         self.Stokes[2],
                         )


        im2 = ax1.pcolormesh(self.Meshes.ThetaVec.Degree,
                             self.Meshes.PhiVec.Degree,
                             self.Stokes[3],
                             shading='auto'
                             )

        cbar = plt.colorbar(im2, ax=ax1, pad=0.15, orientation='horizontal')
        cbar.ax.tick_params(labelsize='small')
        cbar.ax.locator_params(nbins=3)

        ax1.quiver(self.Meshes.ThetaMesh.Degree,
                   self.Meshes.PhiMesh.Degree,
                   self.Stokes[1],
                   self.Stokes[2],
                         )

        plt.show()














# -
