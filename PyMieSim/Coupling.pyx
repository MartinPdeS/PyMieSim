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
sum = np.sum
abs = np.abs

cpdef Coupling(Scatterer, Detector):
     EPhi, EThe = Scatterer.uFarField(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian, 1)
     dOmega     = Detector.Mesh.dOmega.Radian
     Omega      = Detector.Mesh.Omega.Radian
     Filter     = Detector.Filter.Radian
     Scalar     = Detector.Scalar



     if Detector.CouplingMode[1] == 'Centered':
         if Detector.CouplingMode[0] == "Intensity":
             return NoCoherentPointCoupling(Scalar0       = Scalar,
                                            EPhi          = EPhi,
                                            EThe          = EThe,
                                            dOmega        = dOmega,
                                            Filter        = Filter)



         if Detector.CouplingMode[0] == "Amplitude":
             return CoherentPointCoupling(Scalar0       = Scalar,
                                          EPhi          = EPhi,
                                          EThe          = EThe,
                                          dOmega        = dOmega,
                                          Filter        = Filter)



     if Detector.CouplingMode[1] == 'Mean':
         if Detector.CouplingMode[0] == "Intensity":                            # same thing as intensity point coupling
           return NoCoherentPointCoupling(Scalar0       = Scalar,
                                          EPhi          = EPhi,
                                          EThe          = EThe,
                                          dOmega        = dOmega,
                                          Filter        = Filter)



         if Detector.CouplingMode[0] == "Amplitude":
             return CoherentMeanCoupling(Scalar0       = Scalar,
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



cpdef NoCoherentPointCoupling(Scalar0, EPhi, EThe, dOmega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = sum( abs(Scalar0 * EPhi)**2 ) * dOmega * PhiFiltering

    val1 = sum( abs(Scalar0 * EThe)**2 ) * dOmega * TheFiltering

    return val0 + val1



cpdef CoherentPointCoupling(Scalar0, EPhi, EThe, dOmega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = abs( sum( Scalar0 * EPhi ) )**2 * dOmega  * PhiFiltering

    val1 = abs( sum( Scalar0 * EThe ) )**2 * dOmega * TheFiltering

    return val0 + val1



cpdef CoherentMeanCoupling(Scalar0, EPhi, EThe, dOmega, Omega, Filter = None):

    PhiFiltering, TheFiltering = GetFiltering(Filter)

    val0 = sum( abs( Scalar0 * EPhi )**2 * dOmega ) / Omega  * PhiFiltering

    val1 = sum( abs( Scalar0 * EThe )**2 * dOmega ) / Omega  * TheFiltering

    return val0 + val1













#
