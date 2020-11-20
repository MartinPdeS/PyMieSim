#cython: language_level=2
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: initializedcheck=False


from libcpp.vector cimport vector
from libcpp.utility cimport pair
cimport cython

import numpy as np
cimport numpy as np

  

ctypedef double double_t
ctypedef double complex complex128_t


cdef extern from "PyMieCoupling/cpp/MieS1S2.cpp":
    cdef pair[vector[complex128_t], vector[complex128_t]] Cwrapper(double a, double b, vector[double_t] phi);


cpdef tuple MieS1S2(double_t m,
                    double_t x,
                    vector[double_t] phi):

    cdef:
      vector[double_t] n, n2
      vector[complex128_t] S1, S2


    cdef pair[vector[complex128_t], vector[complex128_t]] arr
    arr = Cwrapper(m, x, phi)

    return np.asarray(arr)[0], np.asarray(arr)[1]










# -
