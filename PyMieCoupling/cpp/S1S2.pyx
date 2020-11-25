#cython: language_level=2


from libcpp.vector cimport vector
from libcpp.utility cimport pair
cimport cython
import cython
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF


ctypedef double complex complex128_t


cdef extern from "MieS1S2.cpp":
    cdef pair[vector[complex128_t], vector[complex128_t]] Cwrapper(double a, double b, vector[double] phi);

@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef tuple MieS1S2(double m,
                    double x,
                    phi):


    #cdef pair[vector[complex128_t], vector[complex128_t]]
    arr = Cwrapper(m, x, phi)


    a =  np.asarray(arr)

    return a









# -
