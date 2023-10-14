#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
from numpy import cos, sin, exp, sqrt, pi, linspace, abs, arccos, array, all

from PyMieSim.Physics                 import _Polarization
from PyMieSim.Tools.BaseClasses       import BaseSource
from PyMieSim.Tools.utils             import IO
from PyMieSim.Tools.Constants         import *
from PyMieSim.bin.GaussianBeam        import ( Anm,
                                               Anm_integrand,
                                               Bnm,
                                               Bnm_integrand )


EPS = 1e-20




class PlaneWave(BaseSource):
    """
    .. note::
        Class representing a focuse plane wave beam as a light source for
        light scattering.

    Parameters
    ----------
    Wavelength : float
        Wavelenght of the light field.
    Polarization : float
        Polarization of the light field.
    Amplitude : float
        Maximal value of the electric field at focus point.

    """
    kwargformat = ['Wavelength', 'Polarization']


    def __init__(self,
                 Wavelength   : float,
                 Polarization : float = 0,
                 Amplitude    : float = 1):

        self.GLMT         = False
        self.Wavelength   = Wavelength
        self.k            = 2 * pi / Wavelength
        self.Polarization = Polarization
        self.Amplitude    = Amplitude
        self.H0           = Amplitude / (mu0*c)
        self.H0           = sqrt(eps0/mu0) * self.Amplitude
        self.offset       = asarray([EPS]*3)
        self._BSC_        = None


    def Getidx(self, MaxOrder):
        nlist = []
        mlist = []
        if all( self.offset <= EPS):
            for n in range(*MaxOrder):
                for m in [-1,1]:
                    nlist.append(n)
                    mlist.append(m)

        else:
            for n in range(*MaxOrder):
                for m in range(-n,n+1):
                    nlist.append(n)
                    mlist.append(m)

        return tuple( zip( nlist, mlist ) )


    def GetBSC(self, Precision=None, MaxOrder=30, save=False):

        MaxOrder = (1,MaxOrder)
        idx = self.Getidx(MaxOrder)
        index = pd.MultiIndex.from_tuples(idx, names=["n", "m"])

        BSCTE = r'$BSC_{TE}$'; BSCTM = r'$BSC_{TM}$'
        BSC = pd.DataFrame(columns=[BSCTE, BSCTM], index=index)

        for n, m in idx:

            BSC.at[(n,m), BSCTE] = self.BSC( n, m, mode='TE')
            BSC.at[(n,m), BSCTM] = self.BSC( n, m, mode='TM')

        if save:
            fileName = f"./PyMieSim/BSC/PW_{self.Wavelength}.csv"
            print(f" Saving BSC into file:\n {fileName}")
            BSC.to_csv(f'./{fileName}', mode='w')


        self._BSC_ = BSC.reset_index(level=[0,1]).to_numpy().astype(complex)

        self.MaxOrder = MaxOrder[1]#int(self._BSC_[:,0].max().real)

        return BSC


    def _BSC(self):
        if self._BSC_:
            return self._BSC_
        else:
            print('Need to compute BSC before computing light scattering. [<beam>.GetBSC(Precision=20)]')
            raise




    def BSC(self, n, m, mode='TE'):
        """
        .. note::
            Return the beam shape coefficients
            (:math:`g^{l}_{n, TE}`, :math:`g^{l}_{n, TM}`) for a plane wave.
            (Eq: VI.77 of G&G)

        Parameters
        ----------
        n : :class:`int`
            Order of the expansion.
        m : :class:`int`
            Description of parameter `m`.
        mode : :class:`str`
            Mode of the plane wave, either 'TE' or 'TM'.

        Returns
        -------
        :class:`float`
            Expansion coefficient.

        """
        if m not in [-1,1]: return 0

        if mode == 'TM': return 1 / 2

        if mode == 'TE': return  -m * 1j / 2


    def EField(self, x, y, z):
        """
        .. note::
            Value of the electric field.

            .. math::
                E = E_0 \\exp{\\big( i k (z-z_0) \\big) (p_x \\vec{e_x}, p_y \\vec{e_y})}

        """
        z0 = 0
        Pol = self.Polarization.Radian
        FPhase = exp(1j * self.k * (z - z0))
        return self.Amplitude * FPhase * cos(Pol), \
               self.Amplitude * FPhase * sin(Pol), \
               0

    def HField(self, x, y, z):
        """
        .. note::
            Value of the electric field.

            .. math::
                H = H_0 \\exp{\\big( i k (z-z_0) \\big) (p_x \\vec{e_x}, p_y \\vec{e_y})}

        """
        z0 = 0
        Pol = self.Polarization.Radian
        FPhase = exp(1j * self.k * (z - z0))
        return self.H0 * FPhase * cos(Pol), \
               self.H0 * FPhase * sin(Pol), \
               0


    @property
    def I(self):
        """
        .. note::
            Compute averaged value of Poyinting vector.

            .. math::
                I = <S> = \\frac{c \\epsilon_0}{2} E_0^2

        """

        return c * eps0 / 2 * self.Amplitude**2


    def __repr__(self):

        return IO( f"""
                PlaneWave source beam
                Wavelenght:   {self.Wavelength}
                Polarization: {self.Polarization.Degree}
                """ )


class GaussianBeam(BaseSource):
    """
    .. note::
        Class representing a focuse Gaussian beam as a light source for
        light scattering.

    Parameters
    ----------
    Wavelength : float
        Wavelenght of the light field.
    NA : float
        Numerical aperture of the beam.
    Polarization : float
        Polarization of the light field.
    E0 : float
        Maximal value of the electric field at focus point.
    offset : list
        Distance offset of the beam in comparison to the scatterer.

    """
    def __init__(self,
                 Wavelength   : float,
                 NA           : float,
                 Polarization : float = 0.,
                 E0           : float = 1.,
                 Offset       : list  = [EPS]*3) :

        self.GLMT = True
        Offset = array(Offset); Offset[Offset == 0] = EPS
        self.Wavelength = Wavelength
        self.k = 2 * pi / Wavelength
        self.Polarization = Polarization
        self.E0 = E0
        self.NA = NA
        self.w0 = 2 * self.Wavelength / (pi * NA)
        self.s = 1/(self.k*self.w0)
        self.offset = array(Offset)
        self.Offset = self.offset*self.k
        self.R0 = sqrt(self.Offset[0]**2 + self.Offset[1]**2)
        self.xi = arccos(self.Offset[0]/self.R0)
        self._BSC_ = None



    def Getidx(self, MaxOrder):
        nlist = []
        mlist = []
        if all( self.offset <= EPS):
            for n in range(*MaxOrder):
                for m in [-1,1]:
                    nlist.append(n)
                    mlist.append(m)

        else:
            for n in range(*MaxOrder):
                for m in range(-n,n+1):
                    nlist.append(n)
                    mlist.append(m)

        return tuple( zip( nlist, mlist ) )


    def GetBSC(self, MaxOrder=5, save=False, Sampling=200):
        #MaxOrder = self.GetMaxOrder(Precision)
        MaxOrder = (1,MaxOrder)
        idx = self.Getidx(MaxOrder)
        index = pd.MultiIndex.from_tuples(idx, names=["n", "m"])

        BSCTE = r'$BSC_{TE}$'; BSCTM = r'$BSC_{TM}$'
        BSC = pd.DataFrame(columns=[BSCTE, BSCTM], index=index)

        for n, m in idx:

            BSC.at[(n,m), BSCTE] = self.BSC( n, m, mode='TE')
            BSC.at[(n,m), BSCTM] = self.BSC( n, m, mode='TM')


        if save:
            fileName = f"./PyMieSim/BSC/GB_{self.Wavelength}_{self.NA}.csv"
            print(f" Saving BSC into file:\n {fileName}")
            BSC.to_csv(f'./{fileName}', mode='w')

        self._BSC_ = BSC.reset_index(level=[0,1]).to_numpy().astype(complex)

        self.MaxOrder = MaxOrder[1]#int(self._BSC_[:,0].max().real)

        return BSC

    def _BSC(self):
        if self._BSC_:
            return self._BSC_
        else:
            print('Need to compute BSC before computing light scattering. [<beam>.GetBSC(Precision=20)]')
            raise


    def LoadBSC(dir):
        BSC = pd.read_csv(dir)


    def BSC(self, n, m, mode='TE'):
        """
        .. note::
            Return the beam shape coefficients for a focused Gaussian
            beam using the quadrature method (ref[2]:Eq:17).

            (:math:`g^{l}_{n, TE}`, :math:`g^{l}_{n, TM}`)


        Parameters
        ----------
        n : :class:`int`
            Order of the expansion.
        m : :class:`int`
            Description of parameter `m`.
        mode : :class:`str`
            Mode of the plane wave, either 'TE' or 'TM'.

        Returns
        -------
        :class:`float`
            Beam shape coefficient of order n and degree m.

        """

        if mode == 'TM': return self.Anm(n, m)

        if mode == 'TE': return self.Bnm(n, m)


    def GetMaxOrder(self, Precision=5):
        """
        .. note::
            Method return the range of order to evaluate the BSC's.
            However one should be cautious as a high precision will lead to
            numerical overflow.

        Parameters
        ----------
        Precision : :class:`int`
            Precision parameter (standard is 20)

        Returns
        -------
        :class:`tuple`
            Tuple containing the minimum and maximum order to evaluate.

        """

        termP = sqrt( 2.3 * Precision * ( self.s**(-2) + 4 * self.s**2 * self.Offset[2]**2 ) )
        Orders = (max(1, int(self.R0 - 1/2 - termP) ), int(self.R0 - 1/2 + termP))
        self.Orders = Orders
        self.NOrders = abs(Orders[0]-Orders[1])

        return Orders


    def Phi0(self, x, y, z):
        """
        .. note::
            Method returns:

            .. math::
                \\Phi_0 = -i Q \\exp \\Big[ i Q \\frac{(x-x_0)^2 +
                (y-y_0)^2}{w_0^2} \\Big]
        """

        Q = self.Q(x,y,z)
        term = (x-self.offset[0])**2 + (y-self.offset[1])**2
        term *= 1j * Q / self.w0**2
        term = -1j * Q * exp(term)
        return term



    def Phi0_Spherical(self, r, phi, theta):
        """
        .. note::
            Method returns:

            .. math::
                \\Phi_0 = -i Q \\exp \\Big[ \\frac{i Q}{w_0^2} \\Big(  (r
                \\sin (\\theta) \\cos (\\phi) - x_0^2 ) +
                (r \\sin (\\theta) \\sin (\\phi) - y_0^2 ) \\Big) \\Big]
        """
        Q = self.Q_Spherical(r, theta)
        term0 = -1j * Q
        term1 = 1j * Q / self.w0**2
        term2 = ( r * sin(theta) * cos(phi) - self.offset[0] )**2
        term3 = ( r * sin(theta) * sin(phi) - self.y0 )**2

        return term0 * exp(term1 * ( term2 + term3 ) )


    def Bnm(self, n, m):
        """
        .. note::
            From ref[2]:Eq:18-19
            Method returns:

            .. math::
                B_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}
                (\\beta) - \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_4 \\hat{I}_{m}
                (\\beta) d \\theta

            Where:

            .. math::
                \\gamma_0 &= \\frac{(-i)^n \\rho_n^2}{(2n+1) \\psi_n(\\rho_n)}
                e^{-i Z_0}

                .

                \\gamma_1 &= Q \\exp \\big[ i Q s^2 \\Big( R_0 - \\rho_n
                \\sin (\\theta )  \\big)^2 + i \\rho_n \\cos( \\theta ) \\Big]
                \\hat{P}_n^{|m|} (\\cos ( \\theta ) ) \\sin(\\theta)

                .

                \\gamma_2 &= \\big( 2 Q s^2 \\rho_n \\cos(\\theta) -1 \\big)
                 \\sin(\\theta)

                 .

                \\gamma_3 &= 4 Q s^2 X_0 \\cos(\\theta)

                .

                \\gamma_4 &= 4 Q s^2 Y_0 \\cos(\\theta)

                .

                I_m(\\beta) &= \\frac{1}{2\\pi} \int_0^{2\\pi} e^{z \\cos(\\phi)
                - i m \\phi} d\\phi

                .

                \\tilde{I}_m(\\beta) &= e^{-\\beta - i m \\xi} I_m(\\beta)

                .

                Q_{spherical} &= \\frac{1}{2 ( r \\cos (\\theta) - z_0 )
                /k w_0^2 - i}

        """

        return Bnm(n        = n,
                   m        = m,
                   k        = self.k,
                   w0       = self.w0,
                   Offset   = self.Offset)





    def Anm(self, n, m, Sampling=200):
        """
        .. note::
            From ref[2]:Eq:18-19
            Method returns:

            .. math::
                A_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}
                (\\beta) + \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_3 \\hat{I}_{m}
                (\\beta) d \\theta

            Where:

            .. math::
                \\gamma_0 &= \\frac{(-i)^n \\rho_n^2}{(2n+1)
                \\psi_n(\\rho_n)}e^{-i Z_0}

                .

                \\gamma_1 &= Q \\exp \\big[ i Q s^2 \\Big( R_0 - \\rho_n
                \\sin (\\theta )  \\big)^2 + i \\rho_n \\cos( \\theta ) \\Big]
                \\hat{P}_n^{|m|} (\\cos ( \\theta ) ) \\sin(\\theta)

                .

                \\gamma_2 &= \\big( 2 Q s^2 \\rho_n \\cos(\\theta) -1 \\big)
                \\sin(\\theta)

                .

                \\gamma_3 &= 4 Q s^2 X_0 \\cos(\\theta)

                .

                \\gamma_4 &= 4 Q s^2 Y_0 \\cos(\\theta)

                .

                I_m(\\beta) &= \\frac{1}{2\\pi} \int_0^{2\\pi}
                e^{z \\cos(\\phi) - i m \\phi} d\\phi

                .

                \\tilde{I}_m(\\beta) &= e^{-\\beta - i m \\xi} I_m(\\beta)

                .

                Q_{spherical} &= \\frac{1}{2 ( r \\cos (\\theta) - z_0 )
                /k w_0^2 - i}

        """

        return Anm(n        = n,
                   m        = m,
                   k        = self.k,
                   w0       = self.w0,
                   Offset   = self.Offset)





    def Anm_integrand(self, n, m, Sampling=200):
        """
        .. note::
            Method returns:

            .. math::
                A_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}
                (\\beta) + \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_3 \\hat{I}_{m}
                (\\beta) d \\theta

        """

        return Anm_integrand(n        = n,
                             m        = m,
                             sampling = Sampling,
                             k        = self.k,
                             w0       = self.w0,
                             Offset   = self.Offset)


    def Bnm_integrand(self, n, m, Sampling=200):
        """
        .. note::
            Method returns:

            .. math::
                A_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}
                (\\beta) + \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_3 \\hat{I}_{m}
                (\\beta) d \\theta

        """

        return Bnm_integrand(n        = n,
                             m        = m,
                             sampling = Sampling,
                             k        = self.k,
                             w0       = self.w0,
                             Offset   = self.Offset)




    def EField(self, x, y, z):
        """
        .. note::
            Value of the electric Field (First Order), ref[2]:Eq:6

            .. math::
                E = E_0 \\exp{\\big( i k (z-z_0) \\big)
                (p_x \\vec{e_x}, p_y \\vec{e_y})}
        """

        w0 = 1e-5

        s = 1 / (self.k * self.w0)

        _Q = self.Q(z)

        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * (z - self.offset[2]))

        return self.E0 * self.Phi0(x,y,z) * FPhase * 1, \
               0, \
               self.E0 * self.Phi0(x,y,z) * FPhase * -2*s*_Q*(x-self.offset[0])/self.w0


    def HField(self, x, y, z):
        """
        .. note::
            Value of the electric Field (First Order), ref[2]:Eq:6

            .. math::
                H = H_0 \\exp{\\big( i k (z-z_0) \\big)
                (p_x \\vec{e_x}, p_y \\vec{e_y})}
        """


        _Q = self.Q(x,y,z)
        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * (z - self.offset[2]))

        return 0, \
               self.H0 * self.Phi0(x,y,z) * FPhase * 1, \
               self.H0 * self.Phi0(x,y,z) * FPhase * -2*self.s*_Q*(y-self.y0)/self.w0


    def Er(self, r, phi, theta):
        """
        .. note::
            Value of the electric Field (First Order), ref[2]:Eq:10

            .. math::
                E_r = E_0 \\Phi_0 \\Big[ \\sin (\\theta) \\cos (\\phi)
                - \\frac{2Q}{kw_0^2} \\cos (\\theta) \\big( r \\sin (\\theta)
                \\cos (\\phi) -x_0 \\big)  \\Big] \\exp^{ik (r \\cos
                (\\theta) -z_0 )}
        """

        _Q = self.Q_Spherical(r, theta)

        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * ( r * cos(theta) - self.offset[2] ) )

        term0 = self.E0 * self.Phi0_Spherical(r, phi, theta)

        term1 = sin(theta) * cos(phi)

        term2 =  2 * _Q / ( self.k * self.w0**2 ) * cos(theta)

        term3 = (r * sin(theta) * cos(phi) - self.offset[0])

        Er = term0 * (term1 - term2 * term3) * FPhase

        return Er



    def Hr(self, r, theta, phi):
        """
        .. note::
            Value of the magnetic Field (First Order), ref[2]:Eq:10

            .. math::
                H_r = E_0 \\Phi_0 \\Big[ \\sin (\\theta) \\cos (\\phi) -
                \\frac{2Q}{kw_0^2} \\cos (\\theta) \\big( r \\sin (\\theta)
                \\cos (\\phi) -y_0 \\big)  \\Big] \\exp^{ik (r \\cos
                (\\theta) -z_0 )}
        """

        _Q = self.Q_Spherical(r, theta)

        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * ( r * cos(theta) - self.offset[2] ) )

        term0 = self.H0 * self.Phi0_Spherical(r, phi, theta)

        term1 = sin(theta) * cos(phi)

        term2 =  2 * _Q / ( self.k * self.w0**2 ) * cos(theta)

        term3 = (r * sin(theta) * cos(phi) - self.y0)

        Hr = term0 * (term1 - term2 * term3) * FPhase

        return Hr


    def __repr__(self):

        return IO ( f"""
                Gaussian source beam
                Wavelenght:   {self.Wavelength}
                Polarization: {self.Polarization.Degree}
                NA:           {self.NA}
                Waist (w_0)   {self.w0}
                """ )
