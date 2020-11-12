#cython: language_level=2
cimport cython
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libc.math cimport sqrt, cos, acos, sin, abs, pow


from scipy.special.cython_special cimport jv as jvCython
from scipy.special.cython_special cimport yv as yvCython



ctypedef double double_t
ctypedef double complex complex128_t
ctypedef int int_t

cdef double_t pi = 3.1415926

@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple MieS1S2(double_t m,
                   double_t x,
                   double_t mu):

    cdef int_t nmax = int(2 + x + 4 * pow(x,1/3) )

    n = Arrange(1, nmax+1)

    cdef vector[double_t] n2
    cdef complex128_t S1, S2
    cdef vector[complex128_t] SS1, SS2, an, bn
    cdef int_t i

    for i in range(n.size()):
        n2.push_back( (2 * n[i] + 1) / (n[i] * (n[i] + 1)) )


    if x < 0.5:
        an, bn = LowFrequencyMie_ab(m,x, nmax, n)
    else:
        an, bn = Mie_ab(m,x, nmax, n)

    pin, taun = MiePiTau(mu,nmax)

    cdef int_t lenght = an.size()

    for i in range(lenght):

      SS1.push_back( n2[i] * ( an[i] * pin[i] + bn[i] * taun[i] )  )
      SS2.push_back( n2[i] * (an[i] * taun[i] + bn[i] * pin[i] )  )

    S1 = Sum1(SS1)
    S2 = Sum1(SS2)

    return S1, S2



@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef tuple LowFrequencyMie_ab(double_t m,
                              double_t x,
                              int_t nmax,
                              vector[double_t] n):

    cdef double_t LL, m2, x3, x4, x5, x6
    cdef complex128_t a1, a2, b1, b2
    cdef vector[complex128_t] an, bn
    cdef extern from "complex.h":
        float complex J

    m2 = pow(m, 2)
    LL = ( pow(m,2) - 1 ) / ( pow(m,2) + 2 )
    x3 = pow(x,3.)
    x4 = pow(x,4.)
    x5 = pow(x,5.)
    x6 = pow(x,6.)


    a1 = (-2.*1J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * pow(LL,2.)
    a2 = (-1J * x5 / 15) * (m2 - 1) / (2 * m2 + 3)
    b1 = (-1J * x5 / 45) * (m2 - 1)
    b2 = 0 + 0J

    an.push_back(a1)
    an.push_back(a2)

    bn.push_back(b1)
    bn.push_back(b2)

    return an, bn




@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef tuple Mie_ab(double_t m,
                  double_t x,
                  int_t nmax,
                  vector[double_t] n):

    cdef double_t mx = m*x
    cdef int_t  nmx = int(max(nmax,abs(mx))+16)
    cdef extern from "complex.h":
        float complex J
    cdef vector[complex128_t] gsx, gs1x, an, bn
    cdef vector[double_t] px, chx, p1x, ch1x, Dn, D, da, db


    p1x.push_back( sin(x) )
    ch1x.push_back( cos(x) )

    cdef int_t i
    for i in range(nmax):
      px.push_back(  sqrt(0.5 * pi * x) * jvCython( n[i] + 0.5, x ) )
      chx.push_back(-sqrt(0.5 * pi * x) * yvCython( n[i] + 0.5, x ) )

      p1x.push_back(  sqrt(0.5 * pi * x) * jvCython( n[i] + 0.5, x ) )
      ch1x.push_back(-sqrt(0.5 * pi * x) * yvCython( n[i] + 0.5, x ) )


    for i in range(px.size()):
      gsx.push_back( px[i] - 1J * chx[i] )
      gs1x.push_back( p1x[i] - 1J * ch1x[i] )

    for i in range(nmx):
      Dn.push_back(0)

    for i in range(nmx - 1, 1, -1):
      Dn[i-1] = (i / mx) - (1 / (Dn[i] + i / mx))


    for i in range(nmax):
      D.push_back(Dn[i+1])
      da.push_back( D[i] / m + n[i] / x )
      db.push_back( m * D[i] + n[i] / x )
      an.push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) )
      bn.push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) )

    return an, bn




@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef tuple MiePiTau(double_t mu,
                    int_t nmax):


  cdef vector[double_t] p, t

  p.push_back(1)
  p.push_back(3 * mu)
  t.push_back(mu)
  t.push_back(3.0 * cos(2 * acos(mu) ) )

  cdef unsigned int i

  for i in range(2,int(nmax)):
    p.push_back( ( (2 * i + 1) * ( mu * p[i-1] ) - (i + 1) * p[i-2] ) / i )
    t.push_back( (i + 1) * mu * p[i] - (i + 2) * p[i-1] )

  return p, t











@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef vector[double_t] Arrange(int_t start, int_t end):
    cdef vector[double_t] vec
    cdef int_t i

    for i in range(start, end):
        vec.push_back(i)

    return vec





@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef complex128_t Sum1(vector[complex128_t] arr):
    cdef complex128_t sum = 0.
    cdef int_t size = arr.size()
    cdef int_t i

    for i in range(size):
        sum += arr[i]

    return sum











# -
