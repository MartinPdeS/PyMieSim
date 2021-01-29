
import numpy as np
import fibermodes
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from ai import cs

from PyMieCoupling.classes.BaseClasses import BaseDetector, MeshProperty
from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Source, SMF28, Angle, _Polarization, PlotUnstructuredSphere, interp_at
from PyMieCoupling.physics import FraunhoferDiffraction


class Photodiode(BaseDetector, MeshProperty):
    """Short summary.

    Parameters
    ----------
    NA : float
        Numerical aperture of imaging system.
    Sampling : int
        Number of sampling points for the mode (inside NA).
    GammaOffset : float
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : float
        Angle offset of detector in the direction parallel to polarization.
    Filter : float
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : str
        Methode for computing mode coupling. Either Centered or Mean.
    """
    def __init__(self,
                 NA:                float  = 0.2,
                 Sampling:          int    = 401,
                 GammaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 CouplingMode:      str    = 'Centered'):


        self._CouplingMode = ('Intensity', CouplingMode)

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(NA).Radian,
                                  Sampling    = Sampling,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.GetSpherical()


    def GetSpherical(self):

        self.Scalar = np.ones(self.Meshes.Sampling) #/ np.sqrt((self.Meshes.Sampling))





class LPmode(BaseDetector, MeshProperty):
    """Short summary.

    Parameters
    ----------
    Mode : tuple
        LP mode index l, m.
    NA : float
        Numerical aperture of imaging system.
    Sampling : int
        Number of sampling points for the mode (inside NA).
    InterpSampling : int
        Number of sampling point for interpolation of FarField mode.
    GammaOffset : float
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : float
        Angle offset of detector in the direction parallel to polarization.
    Filter : float
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : str
        Methode for computing mode coupling. Either Centered or Mean.
    """

    def __init__(self,
                 Mode:           tuple,
                 NA:             float = 0.2,
                 Sampling:       int   = 401,
                 InterpSampling: int   = 51,
                 GammaOffset:    float = 0,
                 PhiOffset:      float = 0,
                 Filter:         float =  'None',
                 CouplingMode:   str   = 'Centered',):

        if len(Mode) <= 2: Mode = Mode[0], Mode[1], 'h'

        assert Mode[2] in ['v','h',''], "Mode orientation should either be v [vertical] or h [horizontal]"

        assert CouplingMode in ['Centered','Mean'], "Coupling mode can either be Centered or Mean"

        assert NA < 1, "Numerical aperture has to be under 1 radian"

        self._CouplingMode = ('Amplitude', CouplingMode)

        self._Filter = _Polarization(Filter)

        self.ModeNumber = Mode[0]+1, Mode[1], Mode[2]

        self.InterpSampling = InterpSampling

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(NA).Radian,
                                  Sampling    = Sampling,
                                  PhiOffset   = PhiOffset,
                                  GammaOffset = GammaOffset)

        self.GetFarField()

        self.GetSpherical()



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber[:2]),
                                      940e-9,
                                      CoreDiameter*self.InterpSampling/4,
                                      self.InterpSampling).Ex()

        temp = np.array(temp, copy=False)

        if self.ModeNumber[2] == 'h': temp = temp.T

        self.Cartesian = FraunhoferDiffraction(temp)



    def GetSpherical(self):

        shape = self.Cartesian.shape

        x, y = np.mgrid[-50: 50: complex(shape[0]), -50: 50: complex(shape[1])]

        z = 50 / np.tan(self.Meshes.MaxAngle)

        _, self._phi, self._theta = cs.cart2sp(x.flatten(), y.flatten(), x.flatten()*0+z)




        self.Scalar = interp_at(x           = self._phi.flatten(),
                                y           = self._theta.flatten(),
                                v           = self.Cartesian.astype(np.complex).flatten(),
                                xp          = self.Meshes.base.Phi,
                                yp          = self.Meshes.base.Theta,
                                algorithm   = 'linear',
                                extrapolate = True)




        if False:
            PlotUnstructureData(self.Scalar, self.Meshes.base.Theta, self.Meshes.base.Phi)



    def Plot(self):
        Name = 'Mode Field'
        ThetaMean = np.mean(self.Meshes.Theta.Degree).round(1)
        PhiMean = np.mean(self.Meshes.Phi.Degree).round(1)

        fig, (ax0, ax1) = plt.subplots(1,
                                 2,
                                 figsize=(8,4),
                                 subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea(central_latitude=PhiMean, central_longitude=ThetaMean)})

        im0 = ax0.tricontour(self.Meshes.Theta.Degree,
                             self.Meshes.Phi.Degree,
                             self.Scalar.real,
                             levels=13,
                             linewidths=0.5,
                             colors='k',
                             transform = ccrs.PlateCarree())

        cntr0 = ax0.tricontourf(self.Meshes.Theta.Degree,
                                self.Meshes.Phi.Degree,
                                self.Scalar.real,
                                levels=13,
                                cmap="inferno",
                                transform = ccrs.PlateCarree())


        im1 = ax1.tricontour(self.Meshes.Theta.Degree,
                             self.Meshes.Phi.Degree,
                             self.Scalar.imag,
                             levels=14,
                             linewidths=0.5,
                             colors='k',
                             transform = ccrs.PlateCarree())

        cntr1 = ax1.tricontourf(self.Meshes.Theta.Degree,
                                self.Meshes.Phi.Degree,
                                self.Scalar.imag,
                                levels=14,
                                cmap="inferno",
                                transform = ccrs.PlateCarree())

        plt.colorbar(mappable=cntr1, fraction=0.046, orientation='horizontal', ax=ax1)
        plt.colorbar(mappable=cntr0, fraction=0.046, orientation='horizontal', ax=ax0)


        ax1.plot(self.Meshes.Theta.Degree,
                 self.Meshes.Phi.Degree,
                 'ko',
                 ms=0.1,
                 transform = ccrs.PlateCarree())

        ax0.plot(self.Meshes.Theta.Degree,
                 self.Meshes.Phi.Degree,
                 'ko',
                 ms=0.1,
                 transform = ccrs.PlateCarree())

        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels = False
        gl.bottom_labels = True

        gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False)
        #gl.xlocator = matplotlib.ticker.FixedLocator([])
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels = False
        gl.bottom_labels = True


        ax1.set_title(f'Real Part {Name}')
        ax1.set_ylabel(r'Angle $\phi$ [Degree]')
        ax1.set_xlabel(r'Angle $\theta$ [Degree]')

        ax0.set_title(f'Imaginary Part {Name}')
        ax0.set_ylabel(r'Angle $\phi$ [Degree]')
        ax0.set_xlabel(r'Angle $\theta$ [Degree]')

        #ax1.set_extent([-170, 170, -90, 90], crs=ccrs.PlateCarree())

        #fig.tight_layout()
        plt.show()












# -
