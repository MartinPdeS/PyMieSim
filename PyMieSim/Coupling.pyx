# -*- coding: utf-8 -*-

#cython: language_level=2
#cython: boundscheck=False
#cython: initializedcheck=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: wraparound=False


from libcpp.vector cimport vector
cimport cython
import cython
cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc
from cpython cimport Py_buffer
from libcpp.utility cimport pair


ctypedef double complex complex128_t


cpdef Coupling(Scatterer, Detector):
     FarFieldPara, FarFieldPerp = Scatterer._FarField(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian)

     if Detector.CouplingMode[1] == 'Centered':
         if Detector.CouplingMode[0] == "Intensity":
             Para, Perp = IntensityPointCoupling(Scalar0       = Detector.Scalar,
                                                 Parallel      = FarFieldPara,
                                                 Perpendicular = FarFieldPerp,
                                                 SinMesh       = Detector.Mesh.SinMesh,
                                                 dOmega        = Detector.Mesh.dOmega.Radian,
                                                 Filter        = Detector.Filter.Radian)

             return Para + Perp

         if Detector.CouplingMode[0] == "Amplitude":
             Para, Perp = AmplitudePointCoupling(Scalar0       = Detector.Scalar,
                                                 Parallel      = FarFieldPara,
                                                 Perpendicular = FarFieldPerp,
                                                 SinMesh       = Detector.Mesh.SinMesh,
                                                 dOmega        = Detector.Mesh.dOmega.Radian,
                                                 Filter        = Detector.Filter.Radian)

             return Para + Perp


     if Detector.CouplingMode[1] == 'Mean':

         if Detector.CouplingMode[0] == "Intensity":

             Para, Perp = IntensityMeanCoupling(Scalar0      = Detector.Scalar,
                                               Parallel      = FarFieldPara,
                                               Perpendicular = FarFieldPerp,
                                               SinMesh       = Detector.Mesh.SinMesh,
                                               dOmega        = Detector.Mesh.dOmega.Radian,
                                               Omega         = Detector.Mesh.Omega.Radian,
                                               Filter        = Detector.Filter.Radian)

             return Para + Perp

         if Detector.CouplingMode[0] == "Amplitude":
             Para, Perp = AmplitudeMeanCoupling(Scalar0       = Detector.Scalar,
                                                Parallel      = FarFieldPara,
                                                Perpendicular = FarFieldPerp,
                                                SinMesh       = Detector.Mesh.SinMesh,
                                                dOmega        = Detector.Mesh.dOmega.Radian,
                                                Omega         = Detector.Mesh.Omega.Radian,
                                               Filter        = Detector.Filter.Radian)

             return Para + Perp




cpdef IntensityPointCoupling(Scalar0,
                            Parallel,
                            Perpendicular,
                            SinMesh,
                            dOmega,
                            Filter = None):

    if Filter != None: ParaFiltering = np.cos(Filter)**2; PerpFiltering = np.sin(Filter)**2
    else: ParaFiltering = 1;  PerpFiltering = 1

    Para = np.sum( np.abs(Scalar0 * Parallel)**2 * SinMesh) * dOmega**2

    Perp = np.sum( np.abs(Scalar0 * Perpendicular)**2 * SinMesh) * dOmega**2

    return Para * ParaFiltering, Perp * PerpFiltering


cpdef IntensityMeanCoupling(Scalar0,
                            Parallel,
                            Perpendicular,
                            SinMesh,
                            dOmega,
                            Omega,
                            Filter = None):

    return IntensityPointCoupling(Scalar0       = Scalar0,
                                  Parallel      = Parallel,
                                  Perpendicular = Perpendicular,
                                  SinMesh       = SinMesh,
                                  dOmega        = dOmega,
                                  Filter        = Filter)


cpdef AmplitudePointCoupling(Scalar0,
                             Parallel,
                             Perpendicular,
                             SinMesh,
                             dOmega,
                             Filter = None):

    if Filter != None: ParaFiltering = np.cos(Filter)**2; PerpFiltering = np.sin(Filter)**2
    
    else: ParaFiltering = 1;  PerpFiltering = 1

    Para = np.abs( np.sum( Scalar0 * Parallel * SinMesh) )**2 * dOmega**2

    Perp = np.abs( np.sum( Scalar0 * Perpendicular * SinMesh) )**2 * dOmega**2

    return Para * ParaFiltering, Perp * PerpFiltering



cpdef AmplitudeMeanCoupling(Scalar0,
                            Parallel,
                            Perpendicular,
                            SinMesh,
                            dOmega,
                            Omega,
                            Filter = None):

    if Filter != None: ParaFiltering = np.cos(Filter)**2; PerpFiltering = np.sin(Filter)**2
    else: ParaFiltering = 1;  PerpFiltering = 1

    Para = np.sum( np.abs( Scalar0 * Parallel )**2 * dOmega ) / Omega

    Perp = np.sum( np.abs( Scalar0 * Perpendicular )**2 * dOmega ) / Omega

    return Para * ParaFiltering, Perp * PerpFiltering













#
