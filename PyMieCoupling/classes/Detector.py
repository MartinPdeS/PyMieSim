
import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad, Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Meshes import Meshes
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

class Detector(object):

    def __init__(self,
                 size: float = 50e-6,
                 shape: str = 'circle',
                 wavelength: float = 1e-6,
                 npts: int = 101,
                 ThetaOffset: float = 0,
                 PhiOffset: float = 0):

        self._coupling = 'Intensity'

        self.size, self.wavelength = size, wavelength

        self.ThetaOffset, self.PhiOffset = ThetaOffset, PhiOffset

        self.k = 2 * np.pi / wavelength

        self.npts = npts

        self.GenShift()

        self.GenMeshes()

        self.Field, self.Fourier = self.GenField(shape=shape)


    def GenMeshes(self):

        self.DirectVec = np.linspace(-1.5*self.size/2, 1.5*self.size/2, self.npts)

        self.AngleVec = Direct2Angle(self.DirectVec, self.k)

        self._DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self._AngleBound = np.array( [self.AngleVec[0], self.AngleVec[-1]] )

        self.Meshes = Meshes(npts       = self.npts,
                             ThetaBound = self._AngleBound + self.ThetaOffset,
                             PhiBound   = self._AngleBound + self.PhiOffset)


    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1)*np.pi*np.arange(self.npts)*(self.npts-1)/self.npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        self.shift_grid = shift_grid * shift_grid.T


    def magnificate(self, magnification):

        self.DirectVec /= magnification

        self.GenMeshes()


    def GenField(self, shape: str = 'circle'):

        Field = np.zeros([self.npts, self.npts])

        gridX, gridY = np.meshgrid(self.DirectVec, self.DirectVec)


        if shape == 'square':

            gridX[np.abs(gridX) > self.size/2] = 100

            gridY[np.abs(gridY) > self.size/2] = 100

            grid = gridX + gridY

            Field[grid < 50] = 1

            Field[grid >= 50] = 0


        if shape == 'circle':

            grid = np.sqrt((gridX)**2 + (gridY)**2)

            Field[grid > self.size/2] = 0

            Field[grid <= self.size/2] = 1


        norm = np.sum(np.abs(Field).flatten())

        Field /= norm

        Fourier = Field

        return Field, Fourier


    def GenFigure(self):

        fig = plt.figure(figsize=(15, 5))
        ax0 = fig.add_subplot(131)
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133, projection='polar')

        ax0.set_xlabel(r'Distance X [$\mu$m]')
        ax0.set_ylabel(r'Distance Y [$\mu$m]')
        ax0.set_title('LP mode Field')
        ax0.set_aspect('equal')

        ax1.set_xlabel(r'Angle $\theta$ [deg]')
        ax1.set_ylabel(r'Angle $\phi$ [deg]')
        ax1.set_title('LP mode Far-Field')
        ax1.set_aspect('equal')
        ax1.yaxis.tick_right()

        ax2.set_title('LP polar representation')

        return fig, ax0, ax1, ax2

    def PlotFields(self):

        fig, ax0, ax1, ax2 = self.GenFigure()

        unit = 1e6

        im0 = ax0.pcolormesh(self.DirectVec*unit,
                             self.DirectVec*unit,
                             self.Field,
                             shading='auto')

        im1 = ax1.pcolormesh(self.Meshes.Phi.Vector.Degree,
                             self.Meshes.Theta.Vector.Degree,
                             self.Fourier,
                             shading='auto')

        cbar = fig.colorbar(im0,
                            ax=ax0,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')
        cbar.ax.locator_params(nbins=3)

        cbar = fig.colorbar(im1,
                            ax=ax1,
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        data = (np.abs(self.Fourier)[self.npts//2])

        ax2.plot(self.Meshes.Phi.Vector.Radian, data)

        ax2.fill_between(self.Meshes.Phi.Vector.Radian, min(data), data, color='C0', alpha=0.4)

        plt.show()














# --
