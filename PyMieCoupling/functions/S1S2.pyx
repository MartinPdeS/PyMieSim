import numpy as np
cimport numpy as np
#from scipy.special import jv, yv
from scipy.special.cython_special cimport jv as jvCython
from scipy.special.cython_special cimport yv as yvCython
import cython
from libc.math cimport sqrt, cos, acos, sin, abs


ctypedef double complex complex128_t
ctypedef double double_t
ctypedef long long_t
ctypedef int int_t

global pi
pi = 3.141592


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef MieS1S2(double_t m,
              double_t x,
              double_t mu):

  cdef complex128_t S1, S2
  cdef int_t nmax
  cdef np.ndarray[long_t, ndim=1] n
  cdef np.ndarray[double_t, ndim=1] n2

  nmax = int(2+x+4*(x**(1/3)))
  n = np.arange(1,int(nmax)+1)
  n2 = (2*n+1) / (n*(n+1))

  an, bn = AutoMie_ab(m,x)
  pin, taun = MiePiTau(mu, nmax)

  S1 = ( (n2[0:len(an)] * (an*pin[0:len(an)] + bn * taun[0:len(bn)] ) ) ).sum()

  S2 = (n2[0:len(an)] * (an*taun[0:len(an)] + bn * pin[0:len(bn)] ) ).sum()

  return S1, S2


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef MiePiTau(double_t mu,
              int_t nmax):

  cdef np.ndarray[double_t, ndim=1] p = np.zeros(int(nmax))
  cdef np.ndarray[double_t, ndim=1] t = np.zeros(int(nmax))
  p[0] = 1
  p[1] = 3 * mu
  t[0] = mu
  t[1] = 3.0 * cos(2*acos(mu))

  cdef unsigned int n
  for n in range(2,int(nmax)):
    p[n] = ((2*n+1)*(mu*p[n-1])-(n+1)*p[n-2])/n
    t[n] = (n+1)*mu*p[n]-(n+2)*p[n-1]
  return p, t


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef AutoMieQ(double_t m,
             double_t wavelength,
             double_t diameter,
             complex128_t nMedium = 1.0,
             double_t crossover   = 0.01,
             bint asDict          = False,
             bint asCrossSection  = False):

  cdef double_t nMedium_real = nMedium.real
  cdef double_t m_eff = m / nMedium_real
  cdef double_t wavelength_eff = wavelength / nMedium_real
  cdef double_t x_eff = pi * diameter / wavelength_eff


  if x_eff == 0:
    return 0, 0, 0, 1.5, 0, 0, 0

  elif x_eff < crossover:
    return RayleighMieQ(m, wavelength, diameter, nMedium_real, asDict=asDict, asCrossSection=asCrossSection)

  else:
    return MieQ(m, wavelength, diameter, nMedium_real, asDict=asDict, asCrossSection=asCrossSection)


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef Mie_ab(double_t m,
            double_t x):

  cdef double_t mx = m*x

  cdef int_t nmax = int(2+x+4*(x**(1/3))) #nmax = np.round(2+x+4*(x**(1/3)))

  cdef int_t  nmx = int(max(nmax,abs(mx))+16) #nmx = np.round(max(nmax,np.abs(mx))+16)

  cdef np.ndarray[long_t, ndim=1] n = np.arange(1,nmax+1) #

  cdef np.ndarray[double_t, ndim=1] nu = n + 0.5 #

  cdef double_t sx = sqrt(0.5 * pi * x)

  cdef np.ndarray[complex128_t, ndim=1] px = np.zeros(nmax, dtype='complex')
  
  cdef np.ndarray[complex128_t, ndim=1] chx = np.zeros(nmax, dtype='complex')

  cdef int_t N
  for N in range(nmax):
    px[N] = sx * jvCython(nu[N],x)
    chx[N] = -sx * yvCython(nu[N],x)

  cdef np.ndarray[complex128_t, ndim=1] p1x = np.concatenate([[sin(x)], px[0:nmax-1]]) #p1x = np.append(np.sin(x), px[0:int(nmax)-1])

  cdef np.ndarray[complex128_t, ndim=1] ch1x = np.concatenate([[cos(x)], chx[0:nmax-1]]) #ch1x = np.append(np.cos(x), chx[0:int(nmax)-1])

  cdef np.ndarray[complex128_t, ndim=1] gsx = px-(0+1j)*chx #

  cdef np.ndarray[complex128_t, ndim=1] gs1x = p1x-(0+1j)*ch1x #

  # B&H Equation 4.89
  cdef np.ndarray[complex128_t, ndim=1] Dn = np.zeros(int(nmx),dtype=complex)

  cdef unsigned int i
  for i in range(nmx - 1, 1, -1):
    Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))

  cdef np.ndarray[complex128_t, ndim=1] D = Dn[1:int(nmax)+1] # Dn(mx), drop terms beyond nMax

  cdef np.ndarray[complex128_t, ndim=1] da = D/m+n/x

  cdef np.ndarray[complex128_t, ndim=1] db = m*D+n/x

  cdef np.ndarray[complex128_t, ndim=1] an = (da*px-p1x)/(da*gsx-gs1x)

  cdef np.ndarray[complex128_t, ndim=1] bn = (db*px-p1x)/(db*gsx-gs1x)

  return an, bn


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef RayleighMieQ(double_t m,
                  double_t wavelength,
                  double_t diameter,
                  complex128_t nMedium=1.0,
                  bint asDict=False,
                  bint asCrossSection=False):

  cdef double_t LLabsSq, qratio, g, qsca, qabs, qext, qback, qpr
  cdef complex128_t LL

  cdef double_t nMedium_real = nMedium.real
  m /= nMedium_real
  wavelength /= nMedium_real
  cdef double_t x = pi*diameter/wavelength

  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0

  elif x>0:
    LL = (m**2-1)/(m**2+2) # Lorentz-Lorenz term

    LLabsSq = abs(LL)**2

    qsca = 8*LLabsSq*(x**4)/3 # B&H eq 5.8

    qabs = 4 * x * LL.imag # B&H eq. 5.11

    qext = qsca + qabs

    qback = 1.5 * qsca # B&H eq. 5.9

    qratio = 1.5

    g = 0.

    qpr = qext

    if asCrossSection:
      css = pi*(diameter/2)**2

      cext = css*qext

      csca = css*qsca

      cabs = css*qabs

      cpr = css*qpr

      cback = css*qback

      cratio = css*qratio

      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio




@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def MieQ(double_t     m,
         double_t     wavelength,
         double_t     diameter,
         complex128_t nMedium        = 1.0,
         bint         asDict         = False,
         bint         asCrossSection = False):

  cdef double_t nMedium_real = nMedium.real
  cdef double_t x = pi*diameter/wavelength
  cdef int_t nmax
  cdef double_t x2
  cdef np.ndarray[double_t, ndim=1] n, n1, n2, n3, n4

  wavelength /= nMedium_real

  m /= nMedium_real

  if x == 0:
    return 0, 0, 0, 1.5, 0, 0, 0

  elif x <= 0.05:
    return RayleighMieQ(m, wavelength, diameter, nMedium_real, asDict)

  elif x > 0.05:
    #nmax = np.round(2+x+4*(x**(1/3)))

    nmax = int(2+x+4*(x**(1/3)))

    n = np.arange(1,nmax+1)

    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = Mie_ab(m,x)

    qext = (2/x2) * (n1*(an.real+bn.real)).sum()
    qsca = (2/x2) * (n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2)).sum()
    qabs = qext-qsca

    g1 = [an.real[1:int(nmax)],
          an.imag[1:int(nmax)],
          bn.real[1:int(nmax)],
          bn.imag[1:int(nmax)]]

    g1 = [np.append(x, 0.0) for x in g1]

    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca
    if asCrossSection:
      css = pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def LowFrequencyMieQ(double_t     m,
                     double_t     wavelength,
                     double_t     diameter,
                     complex128_t nMedium        = 1.0,
                     bint         asDict         = False,
                     bint         asCrossSection = False):

  cdef double_t nMedium_real = nMedium.real

  m /= nMedium_real

  wavelength /= nMedium_real

  cdef double_t x = pi*diameter/wavelength

  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    n = np.arange(1,3)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = LowFrequencyMie_ab(m,x)

    qext = (2/x2) * (n1*(an.real+bn.real)).sum()
    qsca = (2/x2) * (n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2)).sum()
    qabs = qext-qsca

    g1 = [an.real[1:2],an.imag[1:2],bn.real[1:2],bn.imag[1:2]]

    g1 = [np.append(x, 0.0) for x in g1]

    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca

    if asCrossSection:
      css = pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio





@cython.boundscheck(False)
@cython.wraparound(False)
def AutoMie_ab(double_t m,
               double_t x):

  if x < 0.5:
    return LowFrequencyMie_ab(m,x)
  else:
    return Mie_ab(m,x)


@cython.boundscheck(False)
@cython.wraparound(False)
def LowFrequencyMie_ab(double_t m,
                       double_t x):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMie_ab
  # B&H page 131
  m2 = m**2
  LL = (m**2-1)/(m**2+2)
  x3 = x**3
  x5 = x**5
  x6 = x**6

  a1 = (-2j*x3/3)*LL-(2j*x5/5)*LL*(m2-2)/(m2+2)+(4*x6/9)*(LL**2)
  a2 = (-1j*x5/15)*(m2-1)/(2*m2+3)
  b1 = (-1j*x5/45)*(m2-1)
  b2 = 0+0j
  an = np.append(a1,a2)
  bn = np.append(b1,b2)
  return an,bn





@cython.nonecheck(True)
@cython.boundscheck(True)
@cython.wraparound(True)
cdef complex128_t CSum(np.ndarray[complex128_t, ndim=1] arr):
  cdef int nrows = arr.shape[0]
  cdef int i
  cdef complex128_t sum
  for i in range(nrows):
    sum += arr[i]

  return sum











  # -
