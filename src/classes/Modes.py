import fibermodes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

from miecoupling.src.functions.converts import rad2deg, deg2rad, Angle2Direct, Direct2Angle


class mode(object):

    def __init__(self, fiber, LPmode, wavelength, npts, magnification=1):

        self.mode = fibermodes.Mode(fibermodes.ModeFamily.HE, *LPmode)

        self.fiber = fiber

        self.wavelength = wavelength

        self.k = 2 * np.pi / wavelength

        self.npts = npts

        self.magnification = magnification

        self.GenShift()

        self.DirectVec = np.linspace(-self.fiber.MaxDirect, self.fiber.MaxDirect, self.npts)

        self.GenMeshes()

        self.Field, self.Fourier = self.GenField()


    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1)*np.pi*np.arange(self.npts)*(self.npts-1)/self.npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        self.shift_grid = shift_grid * shift_grid.T

    def magnificate(self, magnification):

        self.DirectVec /= magnification

        self.GenMeshes()

    def GenMeshes(self):

        self.DirectMesh = np.meshgrid(self.DirectVec, self.DirectVec)

        self.AngleVec = Direct2Angle(self.DirectVec, self.k)

        self.RadVec = deg2rad(self.AngleVec)

        self.MaxAngle = self.AngleVec[0]

        self.AngleMesh = np.meshgrid(self.AngleVec, self.AngleVec)

        self.DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self.AngleBound = [self.AngleVec[0], self.AngleVec[-1]]

    def GenField(self):

        Field = fibermodes.field.Field(self.fiber.source,
                                       self.mode,
                                       fibermodes.Wavelength(self.wavelength),
                                       self.DirectVec[0],
                                       np=self.npts)

        Field = Field.Ex()

        norm = np.sum(np.abs(Field))

        Field /= norm

        Fourier = np.fft.fft2(Field)

        Fourier /= self.shift_grid

        Fourier = np.fft.fftshift(Fourier)

        norm = np.sum(np.abs(Fourier))

        Fourier /= norm

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

        im0 = ax0.pcolormesh(self.DirectVec*unit, self.DirectVec*unit, self.Field)

        im1 = ax1.pcolormesh(self.AngleVec, self.AngleVec, np.imag(self.Fourier))

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

        ax2.plot(self.RadVec, data)

        ax2.fill_between(self.RadVec, min(data), data, color='C0', alpha=0.4)

        plt.show()
