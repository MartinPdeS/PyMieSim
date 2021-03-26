# -*- coding: utf-8 -*-
# http://pymiescatt.readthedocs.io/en/latest/forward.html
import numpy as np
from scipy.special import jv, yv
from scipy.integrate import trapz
import warnings


def Fields(Index,
           Diameter,
           Wavelength,
           nMedium,
           Phi,
           Theta,
           Polarization  = 0,
           E0            = 1,
           R             = 1):

    k = 2*np.pi/Wavelength

    SizeParam = k * (Diameter/2)
    s1s2 = S1S2(Index, Diameter, Wavelength, nMedium, Phi-np.pi/2)

    Parallel = np.outer(s1s2[0], np.sin(Theta) );

    Perpendicular = np.outer(s1s2[1], np.cos(Theta) );

    propagator =  E0 / (k*R) * np.exp(-1j*k*R)

    Parallel *= 1j * propagator

    Perpendicular *= - propagator

    return Parallel, Perpendicular


def coerceDType(d):
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d


def Mie_ab(m,x):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_ab
  mx = m*x
  nmax = np.round(2+x+4*(x**(1/3)))
  nmx = np.round(max(nmax,np.abs(mx))+16)
  n = np.arange(1,nmax+1) #
  nu = n + 0.5 #

  sx = np.sqrt(0.5*np.pi*x)

  px = sx*jv(nu,x) #
  p1x = np.append(np.sin(x), px[0:int(nmax)-1]) #

  chx = -sx*yv(nu,x) #
  ch1x = np.append(np.cos(x), chx[0:int(nmax)-1]) #

  gsx = px-(0+1j)*chx #
  gs1x = p1x-(0+1j)*ch1x #

  # B&H Equation 4.89
  Dn = np.zeros(int(nmx),dtype=complex)
  for i in range(int(nmx)-1,1,-1):
    Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))

  D = Dn[1:int(nmax)+1] # Dn(mx), drop terms beyond nMax
  da = D/m+n/x
  db = m*D+n/x

  an = (da*px-p1x)/(da*gsx-gs1x)
  bn = (db*px-p1x)/(db*gsx-gs1x)


  return an, bn




def LowFrequencyMie_ab(m,x):
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

def AutoMie_ab(m,x):
  if x<0.5:
    return LowFrequencyMie_ab(m,x)
  else:
    return Mie_ab(m,x)





def S1S2(Index, Diameter, Wavelength, nMedium, Phi):
  muList = np.cos(Phi)
  k = 2*np.pi/Wavelength
  SizeParam = k * (Diameter/2)
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieS1S2
  S1, S2 = [], []
  for mu in muList:
      nmax = np.round(2+SizeParam+4*np.power(SizeParam,1/3))
      #print('####',nmax)
      an, bn = AutoMie_ab(Index,SizeParam)
      pin, taun = MiePiTau(mu,nmax)

      n = np.arange(1,int(nmax)+1)
      n2 = (2*n+1)/(n*(n+1))
      S1.append(np.sum(n2[0:len(an)]*(an*pin[0:len(an)]+bn*taun[0:len(bn)])))
      S2.append(np.sum(n2[0:len(an)]*(an*taun[0:len(an)]+bn*pin[0:len(bn)])))

  return S1, S2

def MiePiTau(mu,nmax):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MiePiTau
  p = np.zeros(int(nmax))
  t = np.zeros(int(nmax))
  p[0] = 1
  p[1] = 3*mu
  t[0] = mu
  t[1] = 3.0*np.cos(2*np.arccos(mu))
  for n in range(2,int(nmax)):
    p[n] = ((2*n+1)*(mu*p[n-1])-(n+1)*p[n-2])/n

    t[n] = (n+1)*mu*p[n]-(n+2)*p[n-1]

  return p, t
