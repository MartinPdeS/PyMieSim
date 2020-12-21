
from PyMieCoupling.classes.BaseClasses import BaseFarField
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.classes.Representations import ScalarFarField
import polarTransform
import matplotlib.pyplot as plt
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
        self.Parallel, self.Perpendicular, self.Scalar = (None,)*3
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


        Scalar = polarImageReal + complex(0,1) * polarImageimag

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = Scalar.shape[0],
                                  PhiNpts     = Scalar.shape[1],
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)

        self.Scalar = ScalarFarField(Scalar, Parent=self)




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
        Scalar = np.ones([self.Npts, self.Npts]) / (self.Npts*self.Npts)

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = self.Npts,
                                  PhiNpts     = self.Npts,
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)

        self.Scalar = ScalarFarField(Scalar, Parent=self)






# -
