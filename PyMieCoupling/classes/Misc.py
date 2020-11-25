import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib.ticker as tick
import numpy as cp



global Fontsize, pi
Fontsize, pi = 10, 3.141592

class Source(object):

    def __init__(self,
                 Wavelength: float,
                 Polarization: float):

        self.Wavelength = Wavelength

        self.Polarization = Polarization

        self.k = 2 * pi / Wavelength



class LPField(np.ndarray):

    def __init__(self, array, DirectVec):

        pass


    def __new__(cls, **kwargs):

        this = np.asarray(kwargs['array']).view(cls)

        this.DirectVec = kwargs['DirectVec']

        return this


    def __array_finalize__(self, obj):
        pass


    def __repr__(self):

        txt = 'SHAPE: {0}'.format(self.shape)
        return txt


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
                            self.real,
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
                            self.imag,
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



class LPFourier(np.ndarray):

    def __init__(self, array, Meshes):

        pass



    def __new__(cls, **kwargs):
        this = np.array(kwargs['array'], copy=False)

        this = np.asarray(this).view(cls)

        this.Meshes = kwargs['Meshes']

        return this

    def __array_finalize__(self, obj):
        pass



    def __repr__(self):

        return 'lol'


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
                            self.real,
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
                            self.imag,
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

        ax.set_xlabel('')

        ax.set_ylabel('')

        ax.set_title('Polar representation of far-field collection mode', fontsize = Fontsize)

        ax.plot(self.Meshes.Phi.Vector.Radian, self.__abs__()[self.shape[0]//2])

        ax.fill_between(self.Meshes.Phi.Vector.Radian,
                        0,
                        self.__abs__()[self.shape[0]//2],
                        color='C0',
                        alpha=0.4 )

        plt.show()



# -
