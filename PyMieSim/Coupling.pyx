# -*- coding: utf-8 -*-

#cython: language_level=2
#cython: boundscheck=False
#cython: initializedcheck=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: wraparound=False


cimport cython
import cython
cimport numpy as np
import numpy as np

ctypedef double complex complex128_t


cpdef Coupling(Scatterer, Detector):
     EPhi, EThe = Scatterer.uFarField(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian)
     dOmega     = Detector.Mesh.dOmega.Radian
     Omega      = Detector.Mesh.Omega.Radian
     Filter     = Detector.Filter.Radian
     Scalar     = Detector.Scalar

     if Detector.CouplingMode[1] == 'Centered':
         if Detector.CouplingMode[0] == "Intensity":
             return IntensityPointCoupling(Scalar0       = Scalar,
                                           EPhi          = EPhi,
                                           EThe          = EThe,
                                           dOmega        = dOmega,
                                           Filter        = Filter)



         if Detector.CouplingMode[0] == "Amplitude":
             return AmplitudePointCoupling(Scalar0       = Scalar,
                                           EPhi          = EPhi,
                                           EThe          = EThe,
                                           dOmega        = dOmega,
                                           Filter        = Filter)



     if Detector.CouplingMode[1] == 'Mean':
         if Detector.CouplingMode[0] == "Intensity":                            # same thing as intensity point coupling
           return IntensityPointCoupling(Scalar0       = Scalar,
                                         EPhi          = EPhi,
                                         EThe          = EThe,
                                         dOmega        = dOmega,
                                         Filter        = Filter)



         if Detector.CouplingMode[0] == "Amplitude":
             return AmplitudeMeanCoupling(Scalar0       = Scalar,
                                          EPhi          = EPhi,
                                          EThe          = EThe,
                                          dOmega        = dOmega,
                                          Omega         = Omega,
                                          Filter        = Filter)



cpdef GetFiltering(Filter):

    if Filter != None:
        ParaFiltering = np.cos(Filter)**2
        PerpFiltering = np.sin(Filter)**2
    else:
        ParaFiltering = 1
        PerpFiltering = 1

    return ParaFiltering, PerpFiltering



cpdef IntensityPointCoupling(Scalar0, EPhi, EThe, dOmega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = ( np.sum( np.abs(Scalar0 * EPhi)**2 ) * dOmega ) * PhiFiltering

    val1 = ( np.sum( np.abs(Scalar0 * EThe)**2 ) * dOmega ) * TheFiltering

    return val0 + val1



cpdef AmplitudePointCoupling(Scalar0, EPhi, EThe, dOmega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = ( np.abs( np.sum( Scalar0 * EPhi ) ) * dOmega )**2 * PhiFiltering

    val1 = ( np.abs( np.sum( Scalar0 * EThe ) ) * dOmega )**2 * TheFiltering

    return val0 + val1



cpdef AmplitudeMeanCoupling(Scalar0, EPhi, EThe, dOmega, Omega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = np.sum( np.abs( Scalar0 * EPhi )**2 * dOmega ) / Omega  * PhiFiltering

    val1 = np.sum( np.abs( Scalar0 * EThe )**2 * dOmega ) / Omega  * TheFiltering

    return val0 + val1













#
