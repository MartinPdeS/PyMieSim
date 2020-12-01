#cython: language_level=2


from libcpp.vector cimport vector
cimport cython
import numpy as np
import cython
cimport numpy as np
np.import_array()

ctypedef double complex complex128_t


cdef extern from "MieS1S2.cpp":
    cdef complex128_t*  Cwrapper(double a, double b, vector[double] phi);




@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef np.ndarray MieS1S2(double m,
                         double x,
                         phi):


    cdef complex128_t* array = Cwrapper(m, x, phi)

    cdef np.npy_intp shape[1]

    shape[0] = <np.npy_intp> phi.size * 2

    ndarray = np.PyArray_SimpleNewFromData(1,
                                           shape,
                                           np.NPY_COMPLEX128,
                                           <void *> array
                                            )

    np.PyArray_UpdateFlags(ndarray, ndarray.flags.num | np.NPY_OWNDATA )

    ndarray = np.reshape(ndarray,[2,phi.size])

    return ndarray
