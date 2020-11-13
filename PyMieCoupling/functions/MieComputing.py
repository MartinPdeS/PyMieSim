import numpy as np
from scipy.special import jv, yv



def MieS1S2(m,x,angle):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieS1S2
  #nmax = np.round(2+x+4*np.power(x,1/3))
  MuList = np.cos(angle)
  S1, S2 = [], []
  for mu in MuList:
        nmax = int(2+x+4*(x**(1/3)))
        an, bn = AutoMie_ab(m,x)
        pin, taun = MiePiTau(mu,nmax)
        n = np.arange(1,int(nmax)+1)
        n2 = (2*n+1)/(n*(n+1))
        S1.append( np.sum(n2[0:len(an)]*(an*pin[0:len(an)]+bn*taun[0:len(bn)])) )
        S2.append( np.sum(n2[0:len(an)]*(an*taun[0:len(an)]+bn*pin[0:len(bn)])) )


  return S1, S2


#@jit(nopython=True)
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


#@jit(nopython=True, annotate=True)
def AutoMieQ(m, wavelength, diameter, nMedium=1.0, crossover=0.01, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#AutoMieQ
  nMedium = nMedium.real
  m_eff = m / nMedium
  wavelength_eff = wavelength / nMedium
  x_eff = np.pi*diameter/wavelength_eff
  if x_eff==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x_eff<crossover:
    return RayleighMieQ(m, wavelength, diameter, nMedium, asDict=asDict, asCrossSection=asCrossSection)
  else:
    return MieQ(m, wavelength, diameter, nMedium, asDict=asDict, asCrossSection=asCrossSection)




#@jit(nopython=True)
def RayleighMieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#RayleighMieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    LL = (m**2-1)/(m**2+2) # Lorentz-Lorenz term
    LLabsSq = np.abs(LL)**2
    qsca = 8*LLabsSq*(x**4)/3 # B&H eq 5.8
    qabs=4*x*LL.imag # B&H eq. 5.11
    qext=qsca+qabs
    qback = 1.5*qsca # B&H eq. 5.9
    qratio = 1.5
    g = 0
    qpr = qext
    if asCrossSection:
      css = np.pi*(diameter/2)**2
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





#@jit(nopython=True)
def MieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x<=0.05:
    return RayleighMieQ(m, wavelength, diameter, nMedium, asDict)
  elif x>0.05:
    #nmax = np.round(2+x+4*(x**(1/3)))
    nmax = int(2+x+4*(x**(1/3)))
    n = np.arange(1,nmax+1)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = Mie_ab(m,x)

    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:int(nmax)],
          an.imag[1:int(nmax)],
          bn.real[1:int(nmax)],
          bn.imag[1:int(nmax)]]
    g1 = [np.append(x, 0.0) for x in g1]

    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca
    if asCrossSection:
      css = np.pi*(diameter/2)**2
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


#@jit(nopython=True)
def LowFrequencyMieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    n = np.arange(1,3)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = LowFrequencyMie_ab(m,x)

    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:2],an.imag[1:2],bn.real[1:2],bn.imag[1:2]]

    g1 = [np.append(x, 0.0) for x in g1]

    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca

    if asCrossSection:
      css = np.pi*(diameter/2)**2
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



def Mie_ab(m,x):

  mx = m*x
  #nmax = np.round(2+x+4*(x**(1/3)))
  nmax = int(2+x+4*(x**(1/3)))
  #nmx = np.round(max(nmax,np.abs(mx))+16)
  nmx = int(max(nmax,np.abs(mx))+16)

  n = np.arange(1,nmax+1) #

  nu = n + 0.5 #

  sx = np.sqrt(0.5*np.pi*x)

  px = sx*jv(nu,x) #
  #p1x = np.append(np.sin(x), px[0:int(nmax)-1]) ###############################################________MY____CHANGES

  chx = -sx*yv(nu,x) #
  #ch1x = np.append(np.cos(x), chx[0:int(nmax)-1]) #

  p1x = np.concatenate([[np.sin(x)], px[0:int(nmax)-1]]) #

  ch1x = np.concatenate([[np.cos(x)], chx[0:int(nmax)-1]]) #


  gsx = px-(0+1j)*chx #

  gs1x = p1x-(0+1j)*ch1x #

  # B&H Equation 4.89
  Dn = np.zeros(int(nmx),dtype=complex)

  for i in range(nmx - 1, 1, -1):
    Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))

  D = Dn[1:int(nmax)+1] # Dn(mx), drop terms beyond nMax

  da = D/m+n/x

  db = m*D+n/x

  an = (da*px-p1x)/(da*gsx-gs1x)

  bn = (db*px-p1x)/(db*gsx-gs1x)

  return an, bn



def AutoMie_ab(m,x):
  if x<0.5:
    return LowFrequencyMie_ab(m,x)
  else:
    return Mie_ab(m,x)


#@jit(nopython=True)
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
