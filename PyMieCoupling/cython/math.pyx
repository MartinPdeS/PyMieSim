#cython: language_level=2

cimport cython
cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
from libc.math cimport cos, sin, abs

ctypedef double double_t
ctypedef double complex complex128_t
ctypedef int int_t


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple CCoupling(vector[vector[complex128_t]]     Field,
                      vector[vector[complex128_t]] Parallel,
                      vector[vector[complex128_t]] Perpendicular,
                      vector[vector[double_t]]     Phi,
                      ):

    cdef:
       int_t i, k
       complex128_t SumPerp, SumPara
       double_t dOmega = abs( Phi[0][0] - Phi[0][1] ) * abs( Phi[0][0] - Phi[0][1] )

    for i in range(Field.size()):
        for k in range(Field[0].size()):

            SumPerp += Field[i][k] * Perpendicular[i][k] * abs(sin(Phi[i][k])) * dOmega

            SumPara += Field[i][k] * Parallel[i][k] * abs(sin(Phi[i][k])) * dOmega

    return abs(SumPara) * abs(SumPara) , abs(SumPerp) * abs(SumPerp)



@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef FieldsFromS1S2(np.ndarray[complex128_t, ndim=1] vec1,
                     np.ndarray[double_t, ndim=1] vec2
                     ):
    cdef int_t i, k
    cdef np.ndarray[complex128_t, ndim=2] Parallel = np.zeros((vec1.size, vec2.size), dtype='complex')
    cdef np.ndarray[complex128_t, ndim=2] Perpendicular = np.zeros((vec1.size, vec2.size), dtype='complex')

    for i in range(vec1.size):
      for k in range(vec2.size):
        Parallel[i,k] = vec1[i] * sin(vec2[k])
        Perpendicular[i,k] = vec1[i] * cos(vec2[k])


    return Parallel, Perpendicular


"RfO"



# -
