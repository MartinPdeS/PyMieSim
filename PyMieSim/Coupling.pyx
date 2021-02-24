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
