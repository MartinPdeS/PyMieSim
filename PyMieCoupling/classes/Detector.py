
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import fibermodes

from PyMieCoupling.functions.converts import rad2deg, deg2rad, Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Meshes import Meshes


class DetectorMeta(object):
    """Super() class for Detector and Modes.

    """
    def __init__():
        pass

    def GenShift(self):

        phase_shift = np.exp(-complex(0, 1)*np.pi*np.arange(self.Npts)*(self.Npts-1)/self.Npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        self.shift_grid = shift_grid * shift_grid.T


    def magnificate(self, Magnification):

        self.DirectVec /= Magnification

        self.GenMeshes()


    def GenMeshes(self):

        self.AngleVec = Direct2Angle(self.DirectVec, self.k)

        self._DirectBound = [self.DirectVec[0], self.DirectVec[-1]]

        self._AngleBound = np.array( [self.AngleVec[0], self.AngleVec[-1]] )

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = self._AngleBound + self.ThetaOffset,
                             PhiBound   = self._AngleBound + self.PhiOffset)


class Photodiode(DetectorMeta):
    """Short summary.

    Parameters
    ----------
    size : float
        Size of the detector, [diameter for circle shaped/ side for square].
    shape : str
        Shape of the detector.
    wavelength : float
        Wavelength of the incoming source.
    npts : int
        Number of points defining the rastered meshes.
    ThetaOffset : float
        Offset of theta angle between the detector and source.
    PhiOffset : float
        Offset of phi angle between the detector and source.
    Magnification : float
        Magnification induced by the lense.
    Name : str
        Name of detector [optional for plots].

    Attributes
    ----------
    _name : type
        Description of attribute `_name`.
    _coupling : type
        Description of attribute `_coupling`.
    k : type
        Description of attribute `k`.
    DirectVec : type
        Description of attribute `DirectVec`.
    GenShift : type
        Description of attribute `GenShift`.
    GenMeshes : type
        Description of attribute `GenMeshes`.
    Field : type
        Description of attribute `Field`.
    Fourier : type
        Description of attribute `Fourier`.
    GenField : type
        Description of attribute `GenField`.
    magnificate : type
        Description of attribute `magnificate`.
    size
    wavelength
    ThetaOffset
    PhiOffset
    npts

    """

    def __init__(self,
                 NumericalAperture: float = 0.2,
                 Wavelength: float = 1e-6,
                 Npts: int = 101,
                 ThetaOffset: float = 0,
                 PhiOffset: float = 0,
                 Name: str = 'Detector'):

        self._name = Name

        self._coupling = 'Intensity'

        self.NumericalAperture = NumericalAperture

        self.Wavelength = Wavelength

        self.ThetaOffset, self.PhiOffset = ThetaOffset, PhiOffset

        self.k = 2 * np.pi / Wavelength

        self.Npts = Npts

        self.GenShift()

        self.GenMeshes()

        self.Fourier = self.GenField()


    def GenMeshes(self):
        Angle = np.arcsin(self.NumericalAperture)

        AngleBound = np.array([-Angle, Angle])

        self.Meshes = Meshes(Npts       = self.Npts,
                             ThetaBound = rad2deg(AngleBound) + self.ThetaOffset,
                             PhiBound   = rad2deg(AngleBound) + self.PhiOffset)


    def GenField(self):
        return np.ones( np.shape( self.Meshes.Theta.Mesh.Degree ) )


    def PlotPolar(self):

        fig = plt.figure(figsize=(6,6))

        ax0 = fig.add_subplot(111, projection='polar')

        data = (np.abs(self.Fourier)[self.Npts//2])

        ax0.set_title('Polar representation of far-field collection')

        ax0.plot(self.Meshes.Phi.Vector.Radian, data)

        ax0.fill_between(self.Meshes.Phi.Vector.Radian, 0, data, color='C0', alpha=0.4)

        plt.show()





class LPmode(DetectorMeta):
    """Short summary.

    Parameters
    ----------
    Fiber : float
        The Fiber class used to generate the modes.
    LPmode : str
        The indices of LP mode.
    wavelength : float
        Wavelength of the incoming source.
    npts : int
        Number of points defining the rastered meshes.
    ThetaOffset : float
        Offset of theta angle between the detector and source.
    PhiOffset : float
        Offset of phi angle between the detector and source.
    Magnification : float
        Magnification induced by the lense.
    Name : str
        Name of detector [optional for plots].

    Attributes
    ----------
    _name : type
        Description of attribute `_name`.
    _coupling : type
        Description of attribute `_coupling`.
    k : type
        Description of attribute `k`.
    DirectVec : type
        Description of attribute `DirectVec`.
    GenShift : type
        Description of attribute `GenShift`.
    GenMeshes : type
        Description of attribute `GenMeshes`.
    Field : type
        Description of attribute `Field`.
    Fourier : type
        Description of attribute `Fourier`.
    GenField : type
        Description of attribute `GenField`.
    magnificate : type
        Description of attribute `magnificate`.
    size
    wavelength
    ThetaOffset
    PhiOffset
    npts

    """
    def __init__(self,
                 Fiber: fibermodes.fiber,
                 Mode: tuple,
                 Wavelength: float,
                 Npts: int = 101,
                 Magnification: float = 1.,
                 ThetaOffset: float = 0,
                 PhiOffset: float = 0,
                 Name: str = 'Field detector'):

        self._name, self._coupling, self.Fiber = Name, 'Amplitude', Fiber

        Mode = Mode[0]+1, Mode[1]

        self.Mode = fibermodes.Mode(fibermodes.ModeFamily.HE, *Mode)

        self.Wavelength, self.k = Wavelength, 2 * np.pi / Wavelength

        self.Npts, self.Magnification  = Npts, Magnification

        self.ThetaOffset, self.PhiOffset = ThetaOffset, PhiOffset

        self.DirectVec = np.linspace(-self.Fiber.MaxDirect, self.Fiber.MaxDirect, self.Npts)

        self.GenShift()

        self.GenMeshes()

        self.Field, self.Fourier = self.GenField()

        if Magnification != 1:
            self.magnificate(Magnification)


    def GenField(self):

        Field = fibermodes.field.Field(self.Fiber.source,
                                       self.Mode,
                                       fibermodes.Wavelength(self.Wavelength),
                                       self.DirectVec[0],
                                       np=self.Npts)

        Field = Field.Ex()

        norm = np.sum(np.abs(Field))

        Field /= norm

        Fourier = np.fft.fft2(Field)

        Fourier /= self.shift_grid

        Fourier = np.fft.fftshift(Fourier)

        norm = np.sum(np.abs(Fourier))

        Fourier /= norm

        return Field, Fourier


    def PlotPolar(self):

        fig = plt.figure(figsize=(6,6))

        ax0 = fig.add_subplot(111, projection='polar')

        unit = 1e6

        data = (np.abs(self.Fourier)[self.Npts//2])

        ax0.set_title('[{0}] Polar representation of far-field collection'.format(self._name))

        ax0.plot(self.Meshes.Phi.Vector.Radian, data)

        ax0.fill_between(self.Meshes.Phi.Vector.Radian, 0, data, color='C0', alpha=0.4)

        plt.show()



    def PlotDirectSpace(self):

        unit = 1e6

        fig, axes = plt.subplots(1,2, figsize=(10,6))

        ax0, ax1 = axes[0], axes[1]

        [ ax.set_xlabel(r'Angle $\theta$ [deg]') for ax in axes ]

        [ ax.set_ylabel(r'Angle $\phi$ [deg]') for ax in axes ]

        axes[0].set_title('Real part of LP mode Far-Field')

        ax1.set_title('Imaginary part of LP mode Far-Field')

        im0 = axes[0].pcolormesh(self.DirectVec*unit,
                             self.DirectVec*unit,
                             self.Field,
                             shading='auto')

        im1 = axes[1].pcolormesh(self.Meshes.Phi.Vector.Radian,
                             self.Meshes.Theta.Vector.Radian,
                             np.real(self.Fourier),
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
                            ax=axes[1],
                            orientation="horizontal",
                            pad=0.15,
                            shrink=0.745,
                            format=tick.FormatStrFormatter('%.1e'))

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        plt.show()






# --
