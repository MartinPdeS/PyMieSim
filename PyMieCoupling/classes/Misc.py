import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib.ticker as tick
from PyMieCoupling.functions.converts import deg2rad

global Fontsize
Fontsize = 10

class Source(object):

    def __init__(self,
                 Wavelength: float,
                 Polarization: float):

        self.Wavelength = Wavelength

        self.Polarization = Polarization

        self.k = 2 * np.pi / Wavelength



class LPField(object):
    def __init__(self, input, DirectVec):
        self.Array = input
        self.DirectVec = DirectVec

    def __repr__(self):
        pass


    def Plot(self, Part: str = 'Real'):

        assert Part in ['Real', 'Imag']

        fig, ax = plt.subplots(1,1, figsize=(4,4))

        ax.set_xlabel(r'X Direction [$\mu$m]', fontsize = Fontsize)

        ax.set_ylabel(r'Y Direction [$\mu$m]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')

        if Part == 'Real':
            return self.PlotReal(fig, ax)
        else:
            return self.PlotImag(fig, ax)


    def PlotReal(self, fig, ax):

        ax.set_title('Real part of LP mode Near-Field')

        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            np.real(self.Array),
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax=ax,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        plt.show()


    def PlotImag(self, fig, ax) -> None:

        ax.set_title('Imaginary part of LP mode Near-Field', fontsize = Fontsize)

        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            np.imag(self.Array),
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax=ax,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        plt.show()


class LPFourier(object):
    def __init__(self, input, Meshes):
        self.Array = input
        self.Meshes = Meshes

    def __repr__(self):
        pass


    def Plot(self, Part: str = 'Real') -> None:

        assert Part in ['Real', 'Imag', 'Polar']

        if Part == 'Polar':
            fig, ax = plt.subplots(1,1, figsize=(4,4), subplot_kw={'projection':'polar'})
        else:
            fig, ax = plt.subplots(1,1, figsize=(4,4)
                                   )
        ax.set_xlabel(r'Angle $\theta$ [Degree]', fontsize = Fontsize)

        ax.set_ylabel(r'Angle $\phi$ [Degree]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')

        if Part == 'Real':
            return self.PlotReal(fig, ax)
        elif Part == 'Imag':
            return self.PlotImag(fig, ax)
        elif Part == 'Polar':
            return self.PlotPolar(fig, ax)


    def PlotReal(self, fig, ax) -> None:

        ax.set_title('Real part of LP mode Far-Field', fontsize = Fontsize)

        im0 = ax.pcolormesh(self.Meshes.Phi.Vector.Radian,
                                 self.Meshes.Theta.Vector.Radian,
                                 np.real(self.Array),
                                 shading='auto')

        cbar = fig.colorbar(im0,
                            ax=ax,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        plt.show()



    def PlotImag(self, fig, ax) -> None:

        ax.set_title('Imaginary part of LP mode Far-Field', fontsize = Fontsize)

        im0 = ax.pcolormesh(self.Meshes.Phi.Vector.Radian,
                                 self.Meshes.Theta.Vector.Radian,
                                 np.real(self.Array),
                                 shading='auto')

        cbar = fig.colorbar(im0,
                            ax=ax,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        plt.show()


    def PlotPolar(self, fig, ax):

        data = (np.abs(self.Array)[self.Array.shape[0]//2])

        ax.set_xlabel('')

        ax.set_ylabel('')

        ax.set_title('Polar representation of far-field collection mode', fontsize = Fontsize)

        ax.plot(self.Meshes.Phi.Vector.Radian, data)

        ax.fill_between(self.Meshes.Phi.Vector.Radian, 0, data, color='C0', alpha=0.4)

        plt.show()
