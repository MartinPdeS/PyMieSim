#cython: language_level=2

cimport cython
cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
from libc.math cimport cos, sin

ctypedef double double_t
ctypedef double complex complex128_t
ctypedef int int_t

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








# -
