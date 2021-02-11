#cython: language_level=2
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: initializedcheck=False

cimport cython
from cython.operator cimport dereference as deref
from cython.parallel cimport parallel, prange

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libc.math cimport sqrt, cos, acos, sin, abs, pow, pi


from scipy.special.cython_special cimport jv as jvCython
from scipy.special.cython_special cimport yv as yvCython



ctypedef double double_t
ctypedef double complex complex128_t
ctypedef int int_t


cdef:
    double complex J = 1j



cdef extern from "PyMieCoupling/cpp/MieS1S2.cpp":
    cdef void test(double a, double b, vector[double_t] phi);


cpdef tuple Coupling(double_t                      m,
                     double_t                      x,
                     vector[double_t]&             phi,
                     vector[double_t]&             theta,
                     vector[vector[complex128_t]]& Fourier
                     ):

    cdef:
      vector[vector[complex128_t]] Parallel      = vector[vector[complex128_t]](phi.size(), vector[complex128_t](theta.size()))
      vector[vector[complex128_t]] Perpendicular = vector[vector[complex128_t]](phi.size(), vector[complex128_t](theta.size()))
      vector[complex128_t] S1
      vector[complex128_t] S2
      double_t temp0 = 0, temp1 = 0, *SumParallel = &temp0, *SumPerpendicular = &temp1



    test(m,x, phi)


    MieS1S2(m, x, phi, theta, S1, S2)

    S1S2toFields(S1, S2, theta, Parallel, Perpendicular, SumParallel, SumPerpendicular, Fourier)

    return SumParallel[0], SumPerpendicular[0]





cdef void S1S2toFields(vector[complex128_t]&         S1,
                       vector[complex128_t]&         S2,
                       vector[double_t]&             theta,
                       vector[vector[complex128_t]]& Parallel,
                       vector[vector[complex128_t]]& Perpendicular,
                       double_t*                     SumParallel,
                       double_t*                     SumPerpendicular,
                       vector[vector[complex128_t]]& Fourier,
                       ):

    cdef unsigned long i, k

    for i in range(S1.size()):
        for k in range(theta.size()):
            Parallel[i][k]       = S1[i] * sin(theta[k])
            Perpendicular[i][k]  = S2[i] * cos(theta[k])
            SumParallel[0]       = <double_t>abs( Parallel[i][k] * Fourier[i][k] )
            SumPerpendicular[0]  = <double_t>abs( Perpendicular[i][k] * Fourier[i][k] )




cdef void MieS1S2(double_t               m,
                   double_t              x,
                   vector[double_t]&     phi,
                   vector[double_t]&     theta,
                   vector[complex128_t]& S1,
                   vector[complex128_t]& S2
                   ):

    cdef:
      unsigned long nmax = <unsigned long>(2 + x + 4 * pow(x,1/3) )
      vector[double_t] n, n2
      vector[complex128_t] SS1, SS2, an, bn, pin, taun

    Arrange(1, nmax+1, n, n2)
    print('##############################')

    if x < 0.5:
        LowFrequencyMie_ab(m,x, nmax, n, an, bn)
    else:
        Mie_ab(m,x, nmax, n, an, bn)

    getS1S2(phi, nmax, pin, taun, SS1, SS2, n2, S1, S2, an, bn)






cdef void getS1S2(vector[double_t]&     phi,
                  unsigned long&        nmax,
                  vector[complex128_t]& pin,
                  vector[complex128_t]& taun,
                  vector[complex128_t]& SS1,
                  vector[complex128_t]& SS2,
                  vector[double_t]&     n2,
                  vector[complex128_t]& S2,
                  vector[complex128_t]& S1,
                  vector[complex128_t]& an,
                  vector[complex128_t]& bn
                  ):

  cdef:
    double_t mu
    unsigned long i, k, lenght = <unsigned long>an.size()
    complex128_t temp0 = 0.*J, temp1 = 0.*J;
    complex128_t *SumS1 = &temp0, *SumS2 = &temp1;

  for k in range(phi.size()):

      mu = cos(phi[k])

      MiePiTau(mu, nmax, pin, taun)

      for i in range(lenght):

          SS1.push_back( n2[i] * ( temp0 + MultiplyComplex(bn[i],taun[i] ) ) );



          SS2.push_back( n2[i] * ( MultiplyComplex( an[i], taun[i] ) + MultiplyComplex( bn[i], pin[i] ) ) );


      Sum1(SS1, SumS1)
      Sum1(SS2, SumS2)

      S1.push_back(SumS1[0])
      S2.push_back(SumS2[0])

      SS1.clear()
      SS2.clear()
      pin.clear()
      taun.clear()



cdef void LowFrequencyMie_ab(double_t              m,
                             double_t              x,
                             unsigned long         nmax,
                             vector[double_t]      n,
                             vector[complex128_t]& an,
                             vector[complex128_t]& bn
                             ):

    cdef :
      double_t LL, m2, x3, x4, x5, x6;
      complex128_t a1, a2, b1, b2;

    m2 = m * m;
    LL = ( m2 - 1 ) / ( m2 + 2 );
    x3 = x * x * x;
    x4 = x3 * x;
    x5 = x4 * x;
    x6 = x5 * x;

    a1 = (-2.*1*J * x3 / 3.) * LL - (2.*J * x5 / 5.) * LL * (m2 - 2.) / (m2 + 2.) + (4. * x6 / 9.) * LL * LL;
    a2 = (-1*J * x5 / 15) * (m2 - 1) / (2 * m2 + 3);
    b1 = (-1*J * x5 / 45) * (m2 - 1);
    b2 = 0 + 0*J;

    an.push_back(a1)
    an.push_back(a2)

    bn.push_back(b1)
    bn.push_back(b2)




cdef void Mie_ab(double_t          m,
                 double_t          x,
                 unsigned long     nmax,
                 vector[double_t]& n,
                 vector[complex128_t]& an,
                 vector[complex128_t]& bn
                 ):

    cdef:
      double_t mx = m*x, temp  = sqrt(0.5 * pi * x)
      vector[complex128_t] gsx, gs1x,
      vector[double_t] px, chx, p1x, ch1x, Dn, D, da, db
      unsigned long i, nmx = <unsigned long>(max(nmax,abs(mx))+16)



    p1x.push_back( sin(x) )
    ch1x.push_back( cos(x) )


    for i in range(nmx):
      Dn.push_back(0)

    for i in range(nmx - 1, 1, -1):
      Dn[i-1] = (i / mx) - (1 / (Dn[i] + i / mx))


    for i in range(nmax):
      px.push_back(  temp * jvCython( n[i] + 0.5, x ) )
      chx.push_back(-temp * yvCython( n[i] + 0.5, x ) )

      p1x.push_back(  temp * jvCython( n[i] + 0.5, x ) )
      ch1x.push_back(-temp * yvCython( n[i] + 0.5, x ) )

      gsx.push_back( px[i] - 1*J * chx[i] )
      gs1x.push_back( p1x[i] - 1*J * ch1x[i] )

      D.push_back(Dn[i+1])

      da.push_back( D[i] / m + n[i] / x )
      db.push_back( m * D[i] + n[i] / x )

      an.push_back( (da[i] * px[i] - p1x[i]) / (da[i] * gsx[i] - gs1x[i]) )
      bn.push_back( (db[i] * px[i] - p1x[i]) / (db[i] * gsx[i] - gs1x[i]) )


cdef void MiePiTau(double_t              mu,
                   unsigned long         nmax,
                   vector[complex128_t]& pin,
                   vector[complex128_t]& taun
                   ) nogil:



  pin.push_back(1.)
  pin.push_back(3. * mu)
  taun.push_back(mu)
  taun.push_back(3.0 * cos(2. * acos(mu) ) )

  cdef unsigned long i

  for i in range(2,nmax):
    pin.push_back( ( (2 * i + 1) * ( mu * pin[i-1] ) - (i + 1) * pin[i-2] ) / i )
    taun.push_back( (i + 1) * mu * pin[i] - (i + 2) * pin[i-1] )


cdef void Arrange(unsigned long start,
                  unsigned long end,
                  vector[double_t]& n,
                  vector[double_t]& n2
                  ) nogil:

    cdef unsigned long i

    for i in range(start, end):
        n.push_back(i)
        n2.push_back( (2 * (<double_t>i+1) + 1) / ((<double_t>i+1) * (<double_t>i + 2)) )



cdef void Sum1(vector[complex128_t] arr,
               complex128_t* SumVal
               ) nogil:

    cdef complex128_t val
    SumVal[0] = 0.

    for val in arr:
        SumVal[0] += val



cdef complex128_t MultiplyComplex(complex128_t val0,
                                  complex128_t val1,
                                  ) nogil:

    cdef:
      complex128_t a = val0.real;
      complex128_t b = val0.imag;
      complex128_t c = val1.real;
      complex128_t d = val1.imag;

    c = (a * c + b * d) + J * ( b * c + a * d);

    return c
















# -
