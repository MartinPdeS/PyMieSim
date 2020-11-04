import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib.ticker as tick


global Fontsize, pi
Fontsize, pi = 10, 3.141592

class Source(object):

    def __init__(self,
                 Wavelength: float,
                 Polarization: float):

        self.Wavelength = Wavelength

        self.Polarization = Polarization

        self.k = 2 * pi / Wavelength



class LPField(object):
    def __init__(self, input, DirectVec):
        self.Array = input
        self.DirectVec = DirectVec

    def __repr__(self):
        txt = 'SHAPE: {0}'.format(self.Array.shape)
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

        data, Phi, Theta = CuPy2NumPy(self.Array.real,
                                      self.DirectVec*1e6,
                                      self.DirectVec*1e6)

        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            data,
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

        data, Phi, Theta = CuPy2NumPy(self.Array.imag,
                                      self.DirectVec*1e6,
                                      self.DirectVec*1e6)

        im0 = ax.pcolormesh(self.DirectVec*1e6,
                            self.DirectVec*1e6,
                            data,
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
    def __init__(self,
                 input,
                 Meshes):
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

        data, Phi, Theta = CuPy2NumPy(self.Array.real,
                                      self.Meshes.Phi.Vector.Radian,
                                      self.Meshes.Theta.Vector.Radian)

        im0 = ax.pcolormesh(Phi,
                            Theta,
                            data,
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

        data, Phi, Theta = CuPy2NumPy(self.Array.imag,
                                      self.Meshes.Phi.Vector.Radian,
                                      self.Meshes.Theta.Vector.Radian)

        im0 = ax.pcolormesh(Phi,
                            Theta,
                            data,
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

        data, Phi, Theta = CuPy2NumPy(self.Array.__abs__()[self.Array.shape[0]//2],
                                      self.Meshes.Phi.Vector.Radian,
                                      self.Meshes.Theta.Vector.Radian)

        ax.set_xlabel('')

        ax.set_ylabel('')

        ax.set_title('Polar representation of far-field collection mode', fontsize = Fontsize)

        ax.plot(Phi, data)

        ax.fill_between(Phi, 0, data, color='C0', alpha=0.4)

        plt.show()





class Operation(object):
    def __init__(self):
        pass


    def cos(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.cos(*args)

        elif isinstance(args[0], np.ndarray):
            return np.cos(*args)


    def sin(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.sin(*args)

        elif isinstance(args[0], np.ndarray):
            return np.sin(*args)

    def arcsin(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.arcsin(*args)

        else:
            return np.arcsin(*args)


    def arctan(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.arctan(*args)

        else:
            return np.arctan(*args)        


    def exp(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.exp(*args)

        elif isinstance(args[0], np.ndarray):
            return np.exp(*args)


    def outer(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.outer(args[0], args[1])

        elif isinstance(args[0], np.ndarray):
            return np.outer(args[0], args[1])


    def sqrt(*args):
        if isinstance(args[0], cp.ndarray):
            return cp.sqrt(args)

        elif isinstance(args[0], np.ndarray):
            return np.sqrt(*args)


    def angle(cuda):
        if cuda:
            return cp.angle

        else:
            return np.angle


    def empty(cuda):
        if cuda:
            return cp.empty

        else:
            return np.empty


    def array(cuda):
        if cuda:
            return cp.array

        else:
            return np.array


    def round(cuda):
        if cuda:
            return cp.round

        else:
            return np.round


    def ones(cuda = False):
        if cuda:
            return cp.ones

        else:
            return np.ones


    def fft(cuda = False):
        if cuda:
            return cp.fft

        else:
            return np.fft


    def meshgrid(cuda = False):
        if cuda:
            return cp.meshgrid

        else:
            return np.meshgrid


    def arange(cuda = False):
        if cuda:
            return cp.arange

        else:
            return np.arange


    def linspace(cuda = False):
        if cuda:
            return cp.linspace

        else:
            return np.linspace


    def meshgrid(cuda = False):
        if cuda:
            return cp.meshgrid

        else:
            return np.meshgrid


    pi = 3.141516






def CuPy2NumPy(*items):
    ItemList = []
    for item in items:
        if isinstance(item, cp.ndarray):
            ItemList.append(cp.asnumpy(item))
        else:
            ItemList.append(item)

    temp = tuple(ItemList)

    return (*temp,)




# -
