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
    cdef void* C_GetS1S2(double a,
                         double b,
                         double* phi,
                         int lenght,
                         complex128_t* ptr);

    cdef void* C_GetFields(double a,
                           double b,
                           double* ThetaVec,
                           int Thetalenght,
                           double* PhiVec,
                           int Philenght,
                           complex128_t* Parallel,
                           complex128_t* Perpendicular,
                           double Polarization);


    cdef void* C_GetFieldsNoPolarization(double a,
                                         double b,
                                         double* ThetaVec,
                                         int Thetalenght,
                                         double* PhiVec,
                                         int Philenght,
                                         complex128_t* Parallel,
                                         complex128_t* Perpendicular);



    cdef void* C_GetFieldsFromMesh(double a,
                                  double b,
                                  double* ThetaVec,
                                  double* PhiVec,
                                  int Lenght,
                                  complex128_t* Parallel,
                                  complex128_t* Perpendicular,
                                  double Polarization);


    cdef void* C_GetFieldsNoPolarizationFromMesh(double a,
                                                double b,
                                                double* ThetaVec,
                                                double* PhiVec,
                                                int Lenght,
                                                complex128_t* Parallel,
                                                complex128_t* Perpendicular);


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef GetS1S2(double m,
              double x,
              phi):

    Vector = VectorWrapper(2 * phi.size)
    Vector.add_row()

    cdef:
        np.ndarray[double, ndim=1, mode="c"] phiView = np.asarray(phi, dtype = float, order="C")
        double* phiMesh_ptr = <double *>PyMem_Malloc(sizeof(double*))

    phiMesh_ptr = &phiView[0]

    C_GetS1S2(m,  x, phiMesh_ptr, phi.size, &(Vector.S1S2)[0])

    arr = np.asarray(Vector)

    return np.reshape(arr,[2,phi.size])



@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef GetFields(double m,
                double x,
                ThetaVec,
                PhiVec,
                Polarization):

    cdef:
        np.ndarray[double, ndim=1, mode="c"] ThetaVectorView = np.asarray(ThetaVec, dtype = float, order="C")
        double* ThetaVec_ptr = <double *>PyMem_Malloc(sizeof(double*))

        np.ndarray[double, ndim=1, mode="c"] PhiVectorView = np.asarray(PhiVec, dtype = float, order="C")
        double* PhiVec_ptr = <double *>PyMem_Malloc(sizeof(double*))


    PhiVec_ptr = &PhiVectorView[0]
    ThetaVec_ptr = &ThetaVectorView[0]


    Parallel = VectorWrapper(PhiVec.size * ThetaVec.size)
    Parallel.add_row()

    Perpendicular = VectorWrapper(PhiVec.size * ThetaVec.size)
    Perpendicular.add_row()

    if Polarization is 'None':
      C_GetFieldsNoPolarization(m,
                                x,
                                ThetaVec_ptr,
                                ThetaVec.size,
                                PhiVec_ptr,
                                PhiVec.size,
                                &(Parallel.S1S2)[0],
                                &(Perpendicular.S1S2)[0]);


    else:
      C_GetFields(m,
                  x,
                  ThetaVec_ptr,
                  ThetaVec.size,
                  PhiVec_ptr,
                  PhiVec.size,
                  &(Parallel.S1S2)[0],
                  &(Perpendicular.S1S2)[0],
                  Polarization);


    arr0 = np.asarray(Parallel).reshape([ThetaVec.size, PhiVec.size]).T
    arr1 = np.asarray(Perpendicular).reshape([ThetaVec.size, PhiVec.size]).T

    return arr0, arr1





@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef GetFieldsFromMesh(double m,
                        double x,
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
      C_GetFieldsNoPolarizationFromMesh(m,
                                        x,
                                        ThetaVec_ptr,
                                        PhiVec_ptr,
                                        PhiMesh.size,
                                        &(Parallel.S1S2)[0],
                                        &(Perpendicular.S1S2)[0]);
    else:
      C_GetFieldsFromMesh(m,
                          x,
                          ThetaVec_ptr,
                          PhiVec_ptr, 
                          PhiMesh.size,
                          &(Parallel.S1S2)[0],
                          &(Perpendicular.S1S2)[0],
                          Polarization);


    arr0 = np.asarray(Parallel)
    arr1 = np.asarray(Perpendicular)



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
