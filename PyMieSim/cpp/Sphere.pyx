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

cdef extern from "Mie.cpp":
    cdef void* C_GetS1S2(double        index,
                         double        diameter,
                         double        wavelength,
                         double        nMedium,
                         double*       phi,
                         int           lenght,
                         complex128_t* ptr);

    cdef pair[double, double] C_GetEfficiencies(double a, double b);

    cdef void* Fields(double        index,
                      double        diameter,
                      double        wavelength,
                      double        nMedium,
                      double*       ThetaVec,
                      double*       PhiVec,
                      int           Lenght,
                      complex128_t* Parallel,
                      complex128_t* Perpendicular,
                      double        Polarization);

    cdef void* FieldsNoPolarization(double        index,
                                    double        diameter,
                                    double        wavelength,
                                    double        nMedium,
                                    double*       ThetaVec,
                                    double*       PhiVec,
                                    int           Lenght,
                                    complex128_t* Parallel,
                                    complex128_t* Perpendicular);




cpdef GetS1S2(double index,
              double diameter,
              double wavelength,
              double nMedium,
              phi):


    Vector = VectorWrapper(2 * phi.size)
    Vector.add_row()

    cdef:
        np.ndarray[double, ndim=1, mode="c"] phiView = np.asarray(phi, dtype = float, order="C")
        double* phiMesh_ptr = <double *>PyMem_Malloc(sizeof(double*))

    phiMesh_ptr = &phiView[0]

    C_GetS1S2(index,
              diameter,
              wavelength,
              nMedium,
              phiMesh_ptr,
              phi.size,
              &(Vector.S1S2)[0])

    arr = np.asarray(Vector).reshape([2,phi.size])

    return arr[0], arr[1]


cpdef GetEfficiencies(double m, double x):

    Qsca, Qext = C_GetEfficiencies(m,  x)

    Qabs = Qext - Qsca

    return Qsca, Qext, Qabs



cpdef GetFields(double index,
                double diameter,
                double wavelength,
                double nMedium,
                ThetaMesh,
                PhiMesh,
                Polarization):

    cdef:
        np.ndarray[double, ndim=1, mode="c"] ThetaVectorView = np.asarray(ThetaMesh, dtype = float, order="C")
        double* ThetaVec_ptr = <double *>PyMem_Malloc(sizeof(double*))

        np.ndarray[double, ndim=1, mode="c"] PhiVectorView = np.asarray(PhiMesh, dtype = float, order="C")
        double* PhiVec_ptr = <double *>PyMem_Malloc(sizeof(double*))

    PhiVec_ptr = &PhiVectorView[0]
    ThetaVec_ptr = &ThetaVectorView[0]

    Parallel = VectorWrapper(ThetaMesh.size)
    Parallel.add_row()

    Perpendicular = VectorWrapper(ThetaMesh.size)
    Perpendicular.add_row()

    if Polarization is 'None':
        FieldsNoPolarization(index,
                             diameter,
                             wavelength,
                             nMedium,
                             ThetaVec_ptr,
                             PhiVec_ptr,
                             PhiMesh.size,
                             &(Parallel.S1S2)[0],
                             &(Perpendicular.S1S2)[0]);

    else:
        Fields(index,
               diameter,
               wavelength,
               nMedium,
               ThetaVec_ptr,
               PhiVec_ptr,
               PhiMesh.size,
               &(Parallel.S1S2)[0],
               &(Perpendicular.S1S2)[0],
               Polarization);

    arr0 = np.asarray(Parallel).squeeze()
    arr1 = np.asarray(Perpendicular).squeeze()

    return arr0, arr1





cdef class VectorWrapper:
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
