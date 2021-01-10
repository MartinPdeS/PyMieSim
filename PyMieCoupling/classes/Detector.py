
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from PyMieCoupling.classes.BaseClasses import BaseDetector
import fibermodes
from scipy.interpolate import griddata
import polarTransform

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.classes.Representations import ScalarFarField
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Source, SMF28, Angle, _Polarization, PlotUnstructureData



global cmap
cmap = 'RdBu'


class Photodiode(BaseDetector):

    def __init__(self,
                 Source:            Source = None,
                 NA:                float  = 0.2,
                 Samples:           int    = 401,
                 GammaOffset:       float  = 0,
                 PhiOffset:         float  = 0,
                 Filter:            float  = 'None',
                 Name:              str    = 'Intensity Detector'):

        self._name, self._coupling  = Name, 'Intensity'

        self._NA = NA

        self.Source, self.Samples = Source, Samples

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.Meshes = AngleMeshes(MaxAngle    = NA2Angle(self._NA).Radian,
                                  Samples     = self.Samples,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)

        self.GetSpherical()


    def GetSpherical(self):

        self.Samples = self.Meshes.Phi.Radian.shape

        Scalar = np.ones(self.Samples) / (self.Samples)

        self.Scalar = ScalarFarField(Scalar, Parent=self)


class LPmode(BaseDetector):

    def __init__(self,
                 Mode:          tuple,
                 Source:        Source,
                 Orientation:   str   = 'h',
                 NA:            float = 0.2,
                 Samples:       int   = 401,
                 GammaOffset:   float = 0,
                 PhiOffset:     float = 0,
                 Filter:        float =  'None',
                 Name:          str   = 'Amplitude detector'):

        assert Orientation in ['v','h'], "Orientation should either be v [vertical] or h [horizontal]'"

        self._name, self._coupling = Name, 'Amplitude'

        self._NA = NA

        self.MaxAngle = NA2Angle(self._NA).Radian

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.ModeNumber = Mode[0]+1, Mode[1]

        self.Source, self.Samples, self.Npts = Source, Samples, int(Samples/8)*2 + 1

        self.Orientation = Orientation

        self.debug = False

        MaxAngle = NA2Angle(self._NA).Radian

        self.Meshes = AngleMeshes(MaxAngle    = MaxAngle,
                                  Samples     = self.Samples,
                                  PhiOffset   = self._PhiOffset,
                                  GammaOffset = self._GammaOffset)


        self.GetFarField()

        self.GetSpherical()



    def GetFarField(self):

        Fiber, CoreDiameter = SMF28()

        temp = fibermodes.field.Field(Fiber.source,
                                      fibermodes.Mode(fibermodes.ModeFamily.HE, *self.ModeNumber),
                                      self.Source.Wavelength,
                                      CoreDiameter*self.Npts/6,
                                      self.Npts).Ex()

        temp = np.array(temp, copy=False)

        #temp /= (temp.__abs__()).sum()

        if self.Orientation == 'h': temp = temp.T

        temp = np.fft.fft2(temp)

        temp /= self.GenShift()

        self.Cartesian  = np.fft.fftshift(temp)


    def GenShift(self):
        phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(self.Npts)*(self.Npts-1)/self.Npts)

        shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

        return shift_grid * shift_grid.T




    def GetSpherical(self):
        polarImageReal, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.real,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = self.Cartesian.real.shape[0]//2,
                                                                        finalAngle    = 2*np.pi)

        polarImageimag, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.imag,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = self.Cartesian.real.shape[0]//2,
                                                                        finalAngle    = 2*np.pi)

        self.Polar = polarImageReal + complex(0,1) * polarImageimag

        shape = self.Polar.imag.shape

        print('####', np.pi/2-self.MaxAngle, np.max(self.Meshes.base.Phi.__abs__()))

        ThetaMesh, PhiMesh = np.mgrid[-np.pi: np.pi:complex(shape[0]),
                                      np.pi/2: np.pi/2-self.MaxAngle*1.01 :complex(shape[1])]


        self.Samples = self.Meshes.Phi.Radian.shape

        Scalar = griddata((PhiMesh.flatten(), ThetaMesh.flatten()),
                           self.Polar.astype(np.complex).flatten(),
                           (self.Meshes.base.Phi, self.Meshes.base.Theta),
                           #fill_value = np.nan,
                           method     = 'nearest')

        norm = np.sqrt( np.sum((self.Meshes.SinMesh * Scalar.__abs__())**2) )

        Scalar /=  norm

        if self.debug:
            PlotUnstructureData(Scalar, self.Meshes.base.Theta, self.Meshes.base.Phi)

        self.Scalar = ScalarFarField(Scalar, Parent=self)



    def Plot(self, num=600):

        phi, theta = np.mgrid[-np.pi/2:np.pi/2:complex(0,num),
                              -np.pi:np.pi:complex(0,num)]

        ziReal = griddata((self.Meshes.Theta.Radian,
                           self.Meshes.Phi.Radian),
                           self.Scalar.astype(np.complex),
                           (theta.flatten(), phi.flatten()),
                           fill_value = complex(np.nan, np.nan),
                           method     = 'linear')


        fig, ax = plt.subplots(nrows      = 1,
                               ncols      = 2,
                               figsize    = (8, 3),
                               subplot_kw = {'projection':'mollweide'}
                               )

        im0 = ax[0].pcolormesh(theta,
                               phi,
                               ziReal.real.reshape([num,num]),
                               cmap=cmap,
                               shading='auto')

        cbar = plt.colorbar(mappable=im0, orientation='horizontal', ax=ax[0], format='%.0e')
        cbar.ax.tick_params(labelsize='small')

        ax[0].set_title('Real Part\n Far-Field Spherical Coordinates')
        ax[0].set_ylabel(r'Angle $\phi$ [Degree]')
        ax[0].set_xlabel(r'Angle $\theta$ [Degree]')
        ax[0].grid()

        im1 = ax[1].pcolormesh(theta,
                               phi,
                               ziReal.imag.reshape([num,num]),
                               cmap=cmap,
                               shading='auto')
        cbar = plt.colorbar(mappable=im1, orientation='horizontal', ax=ax[1], format='%.0e')
        cbar.ax.tick_params(labelsize='small')

        ax[1].set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax[1].set_xlabel(r'Angle $\theta$ [Degree]')
        ax[1].grid()

        ax[0].scatter(self.Meshes.Theta.Radian, self.Meshes.Phi.Radian, s=0.1,c='k')
        ax[1].scatter(self.Meshes.Theta.Radian, self.Meshes.Phi.Radian, s=0.1,c='k')

        fig.tight_layout()





























# --
