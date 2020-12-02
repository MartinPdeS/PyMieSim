#cython: language_level=2


from libcpp.vector cimport vector
cimport cython
import cython
cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc
from cpython cimport Py_buffer


ctypedef double complex complex128_t

cdef extern from "MieS1S2.cpp":
    cdef void* Vec_Cwrapper(double a,
                            double b,
                            double* phi,
                            int lenght,
                            complex128_t* ptr);


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef GetS1S2(double m,
              double x,
              phi):

    M = Matrix(2 * phi.size)

    M.add_row()

    cdef:
        np.ndarray[double, ndim=1, mode="c"] a_cython = np.asarray(phi, dtype = float, order="C")
        double* point_to_a = <double *>PyMem_Malloc(sizeof(double*))

    point_to_a = &a_cython[0]

    Vec_Cwrapper(m,  x, point_to_a, phi.size, &(M.S1S2)[0])

    arr = np.asarray(M)

    return np.reshape(arr,[2,phi.size])



cdef class Matrix:
    cdef:
        Py_ssize_t ncols
        Py_ssize_t shape[2]
        Py_ssize_t strides[2]
        vector[complex128_t] S1S2
        int view_count

    def __cinit__(self, Py_ssize_t ncols):
        self.ncols = ncols
        self.view_count = 0

    def add_row(self):
        """Adds a row, initially zero-filled."""
        if self.view_count > 0:
            raise ValueError("can't add row while being viewed")

        self.S1S2.resize(self.S1S2.size() + self.ncols)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.S1S2[0])

        self.shape[0] = self.S1S2.size() / self.ncols
        self.shape[1] = self.ncols

        self.strides[1] = <Py_ssize_t>(  <char *>&(self.S1S2[1]) - <char *>&(self.S1S2[0]))

        self.strides[0] = self.ncols * self.strides[1]

        buffer.buf = <char *>&(self.S1S2[0])
        buffer.format = 'Zd'                     # float
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.S1S2.size() * itemsize
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1


# -
