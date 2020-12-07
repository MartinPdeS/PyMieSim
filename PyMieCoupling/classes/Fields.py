import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib.ticker as tick
import numpy as cp



global Fontsize, pi, cmapPad
Fontsize, pi, cmapPad = 7, 3.141592, 0.2

class Source(object):

    def __init__(self,
                 Wavelength:   float,
                 Polarization: float,
                 Power:        float = 1):

        self.Wavelength = Wavelength
        if Polarization:
            self.Polarization = Angle(Polarization)
        else:
            self.Polarization = None

        self.k = 2 * pi / Wavelength

        self.Power = Power



class LPField(np.ndarray):

    def __init__(self, array, DirectVec):

        pass


    def __new__(cls, **kwargs):

        this = np.asarray(kwargs['array']).view(cls)

        this.DirectVec = kwargs['DirectVec']

        return this


    def __array_finalize__(self, obj):
        pass


    def Plot(self):

        fig, ax = plt.subplots(1,2, figsize=(6,3))

        self.PlotReal(fig, ax[0])

        self.PlotImag(fig, ax[1])

        plt.show()


    def PlotReal(self, fig, ax):

        ax.set_title('Real part of LP mode Near-Field', fontsize = Fontsize)

        ax.set_xlabel(r'X Direction [$\mu$m]', fontsize = Fontsize)

        ax.set_ylabel(r'Y Direction [$\mu$m]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')


        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            self.real,
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax          = ax,
                            orientation = "horizontal",
                            pad         = cmapPad,
                            shrink      = 0.745,
                            format      = tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize=7)

        cbar.ax.locator_params(nbins=3)


    def PlotImag(self, fig, ax) -> None:

        ax.set_title('Imaginary part of LP mode Near-Field', fontsize = Fontsize)

        ax.set_xlabel(r'X Direction [$\mu$m]', fontsize = Fontsize)

        ax.set_ylabel(r'Y Direction [$\mu$m]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')


        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            self.imag,
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax          = ax,
                            orientation = "horizontal",
                            pad         = cmapPad,
                            shrink      = 0.745,
                            format      = tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize=7)

        cbar.ax.locator_params(nbins=3)




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




    def Plot(self, Part: str = 'Real') -> None:

        if Part == 'Polar':
            fig, ax = plt.subplots(1,1, figsize=(3,3), subplot_kw={'projection':'polar'})
            self.PlotPolar(fig, ax)

        else:
            fig, ax = plt.subplots(1,2, figsize=(6,3))
            self.PlotReal(fig, ax[0])

            self.PlotImag(fig, ax[1])

        plt.show()



    def PlotReal(self, fig, ax) -> None:

        ax.set_xlabel(r'Angle $\theta$ [Degree]', fontsize = Fontsize)

        ax.set_ylabel(r'Angle $\phi$ [Degree]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')

        ax.set_title('Real part of LP mode Far-Field', fontsize = Fontsize)

        im0 = ax.pcolormesh(self.Meshes.Phi.Vector.Degree,
                            self.Meshes.Theta.Vector.Degree,
                            self.real,
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax          = ax,
                            orientation = "horizontal",
                            pad         = cmapPad,
                            shrink      = 0.745,
                            format      = tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize=7)

        cbar.ax.locator_params(nbins=3)




    def PlotImag(self, fig, ax) -> None:

        ax.set_xlabel(r'Angle $\theta$ [Degree]', fontsize = Fontsize)

        ax.set_ylabel(r'Angle $\phi$ [Degree]', fontsize = Fontsize)

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')

        ax.set_title('Imaginary part of LP mode Far-Field', fontsize = Fontsize)

        im0 = ax.pcolormesh(self.Meshes.Phi.Vector.Degree,
                            self.Meshes.Theta.Vector.Degree,
                            self.imag,
                            shading='auto')

        cbar = fig.colorbar(im0,
                            ax          = ax,
                            orientation = "horizontal",
                            pad         = cmapPad,
                            shrink      = 0.745,
                            format      = tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize=7)

        cbar.ax.locator_params(nbins=3)



    def PlotPolar(self, fig, ax):

        ax.tick_params(labelsize='small')

        ax.set_aspect('equal')

        ax.set_xlabel('')

        ax.set_ylabel('')

        ax.set_title('Polar representation of far-field collection mode', fontsize = Fontsize)

        ax.plot(self.Meshes.Phi.Vector.Radian, self.__abs__()[self.shape[0]//2])

        ax.fill_between(self.Meshes.Phi.Vector.Radian,
                        0,
                        self.__abs__()[self.shape[0]//2],
                        color='C0',
                        alpha=0.4 )






class Angle(object):

    def __init__(self, input):
        self.Degree = input
        self.Radian = deg2rad(input)



# -
