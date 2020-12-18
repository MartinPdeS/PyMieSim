
from PyMieCoupling.classes.BaseClasses import BaseFarField
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.functions.converts import NA2Angle

import polarTransform
import numpy as np

class LPFarField(BaseFarField):

    def __init__(self,
                 Input,
                 Size:        float,
                 Npts:        int     = 101,
                 NA:          float   = 0.2,
                 PhiOffset:   float   = 0,
                 ThetaOffset: float   = 0):

        self.Cartesian = Input
        self.Size, self.Npts, self._NA  = Size, Npts, NA
        self._PhiBound, self._ThetaBound  =  np.asarray( [0, NA2Angle(self._NA)] ), np.asarray([-180, 180])
        self._PhiOffset, self._ThetaOffset = PhiOffset, ThetaOffset
        self.GetSpherical()

    def GetSpherical(self):
        polarImageReal, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.real,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        polarImageimag, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.imag,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        self.Spherical = polarImageReal + complex(0,1) * polarImageimag

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = self.Spherical.shape[0],
                                  PhiNpts     = self.Spherical.shape[1],
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)


class LPNearField(object):

    def __init__(self, Input, Size, Npts):
        self.Cartesian = Input
        self.Size = Size
        self.Npts = Npts

        self.Meshes = DirectMeshes(Npts   = self.Npts,
                                  XBound = [-self.Size/2, self.Size/2],
                                  YBound = [-self.Size/2, self.Size/2],
                                  XNpts  = self.Npts,
                                  YNpts  = self.Npts)


    def Plot(self):
        fig = plt.figure(figsize=(6,3))
        ax0 = fig.add_subplot(121)
        ax1 = fig.add_subplot(122)

        ax0.pcolormesh(self.Meshes.X.Vector,
                       self.Meshes.Y.Vector,
                       self.Cartesian.real,
                       shading='auto')

        ax0.set_title('Real Part \n Near-Field Cartesian Coordinates')
        ax0.set_xlabel(r'X-Distance x [$\mu$m]')
        ax0.set_ylabel('Y-Distance y  [$\mu$m]')

        ax1.pcolormesh(self.Meshes.X.Vector,
                       self.Meshes.Y.Vector,
                       self.Cartesian.imag,
                       shading='auto')

        ax1.set_title('Imaginary Part \n Near-Field Cartesian Coordinates')
        ax1.set_xlabel(r'X-Distance x  [$\mu$m]')
        ax1.set_ylabel(r'Y-Distance y  [$\mu$m]')
        fig.tight_layout()


class Detector_FarField(BaseFarField):

    def __init__(self,
                 Npts:        int   = 101,
                 NA:          float = 0.2,
                 ThetaOffset: float = 0,
                 PhiOffset:   float = 0):


        self.Npts, self._NA = Npts, NA

        self._PhiBound, self._ThetaBound  =  np.asarray( [0, NA2Angle(self._NA)] ), np.asarray([-180, 180])

        self._PhiOffset, self._ThetaOffset = PhiOffset, ThetaOffset

        self.GetSpherical()



    def GetSpherical(self):
        self.Spherical = np.ones([self.Npts, self.Npts]) / (self.Npts*self.Npts)

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = self.Npts,
                                  PhiNpts     = self.Npts,
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)






# -
