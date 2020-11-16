#cython: language_level=2
cimport cython
from cython.operator cimport dereference as deref
from cython.parallel cimport parallel, prange

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libc.math cimport sqrt, cos, acos, sin, abs, pow


from scipy.special.cython_special cimport jv as jvCython
from scipy.special.cython_special cimport yv as yvCython



ctypedef double double_t
ctypedef double complex complex128_t
ctypedef int int_t
ctypedef vector[complex128_t] CVec
ctypedef vector[CVec] CMatrix


cdef double_t pi = 3.1415926

@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple MieS1S2(double_t m,
                    double_t x,
                    vector[double_t] phi,
                    vector[double_t] theta):

    cdef extern from "complex.h":
        float complex J

    cdef:
      Py_ssize_t nmax = int(2 + x + 4 * pow(x,1/3) )
      vector[double_t] n, n2
      vector[complex128_t] SS1, SS2, an, bn, pin, taun, S1, S2
      Py_ssize_t i, k
      complex128_t temp = 0.*J, *SumS1 = &temp, *SumS2 = &temp


    Arrange(1, nmax+1, n)

    for i in range(n.size()):
        n2.push_back( (2 * n[i] + 1) / (n[i] * (n[i] + 1)) )


    if x < 0.5:
        LowFrequencyMie_ab(m,x, nmax, n, an, bn)
    else:
        Mie_ab(m,x, nmax, n, an, bn)

    cdef:
        unsigned int lenght = an.size()
        double_t mu


    for k in range(phi.size()):

        mu = cos(phi[k])

        MiePiTau(mu,nmax, pin, taun)

        for i in range(lenght):

            SS1.push_back( n2[i] * ( an[i] * pin[i] + bn[i] * taun[i] )  )
            SS2.push_back( n2[i] * (an[i] * taun[i] + bn[i] * pin[i] )  )

        Sum1(SS1, SumS1)
        Sum1(SS2, SumS2)


        S1.push_back(SumS1[0])
        S2.push_back(SumS2[0])

        SS1.clear()
        SS2.clear()
        pin.clear()
        taun.clear()

    return S1, S2



@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void LowFrequencyMie_ab(double_t m,
                             double_t x,
                             Py_ssize_t nmax,
                             vector[double_t] n,
                             vector[complex128_t]& an,
                             vector[complex128_t]& bn
                             ):

    cdef :
      double_t LL, m2, x3, x4, x5, x6
      complex128_t a1, a2, b1, b2
      extern from "complex.h":
          float complex J

    m2 = m*m
    LL = ( m2 - 1 ) / ( m2 + 2 )
    x3 = x*x*x
    x4 = x3*x
    x5 = x4*x
    x6 = x5*x


    a1 = (-2.*1J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL*LL
    a2 = (-1J * x5 / 15) * (m2 - 1) / (2 * m2 + 3)
    b1 = (-1J * x5 / 45) * (m2 - 1)
    b2 = 0 + 0J


    an.push_back(a1)
    an.push_back(a2)

    bn.push_back(b1)
    bn.push_back(b2)




@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Mie_ab(double_t m,
                 double_t x,
                 Py_ssize_t nmax,
                 vector[double_t] n,
                 vector[complex128_t]& an,
                 vector[complex128_t]& bn
                 ):

    cdef:
      double_t mx = m*x,
      double_t temp  = sqrt(0.5 * pi * x)
      int_t  nmx = int(max(nmax,abs(mx))+16)
      vector[complex128_t] gsx, gs1x,
      vector[double_t] px, chx, p1x, ch1x, Dn, D, da, db
      Py_ssize_t i
      extern from "complex.h":
        float complex J


    p1x.push_back( sin(x) )
    ch1x.push_back( cos(x) )


    for i in range(nmx):
      Dn.push_back(0)

    for i in range(nmx - 1, 1, -1):
      Dn[i-1] = (i / mx) - (1 / (Dn[i] + i / mx))


    with nogil:
      for i in range(nmax):
        px.push_back(  temp * jvCython( n[i] + 0.5, x ) )
        chx.push_back(-temp * yvCython( n[i] + 0.5, x ) )

        p1x.push_back(  temp * jvCython( n[i] + 0.5, x ) )
        ch1x.push_back(-temp * yvCython( n[i] + 0.5, x ) )

        gsx.push_back( px[i] - 1J * chx[i] )
        gs1x.push_back( p1x[i] - 1J * ch1x[i] )

        D.push_back(Dn[i+1])

        da.push_back( D[i] / m + n[i] / x )
        db.push_back( m * D[i] + n[i] / x )

        an.push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) )
        bn.push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) )





@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void MiePiTau(double_t mu,
                   int_t nmax,
                   vector[complex128_t]& pin,
                   vector[complex128_t]& taun
                   ) nogil:



  pin.push_back(1.)
  pin.push_back(3. * mu)
  taun.push_back(mu)
  taun.push_back(3.0 * cos(2. * acos(mu) ) )

  cdef Py_ssize_t i

  for i in prange(2,nmax, num_threads=1):
    pin.push_back( ( (2 * i + 1) * ( mu * pin[i-1] ) - (i + 1) * pin[i-2] ) / i )
    taun.push_back( (i + 1) * mu * pin[i] - (i + 2) * pin[i-1] )











@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Arrange(Py_ssize_t start,
                  Py_ssize_t end,
                  vector[double_t]& n) nogil:
    cdef Py_ssize_t i

    for i in range(start, end):
        n.push_back(i)






@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Sum1(vector[complex128_t] arr,
               complex128_t* SumVal) nogil:

    cdef complex128_t val
    SumVal[0] = 0.

    for val in arr:
        SumVal[0] += val























# -
