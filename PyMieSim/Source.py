#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.integrate import trapz as integrate
from scipy.special import spherical_jn as jn
from numpy import cos, sin, exp, sqrt, pi, linspace, abs, arccos, array

from PyMieSim.Physics import _Polarization
from PyMieSim.BaseClasses import BaseSource
from PyMieSim.Constants import eps0, mu0
from PyMieSim.Special import Xi, Pnm_, r8_factorial

class PlaneWave(BaseSource):
    def __init__(self,
                 Wavelength:   float,
                 Polarization: float = 0,
                 E0:           float = 1):

        self.Wavelength = Wavelength
        self.k = 2 * pi / Wavelength
        self.Polarization = _Polarization(Polarization)
        self.E0 = E0
        self.H0 = sqrt(eps0/mu0) * self.E0
        self.x0, self.y0, self.z0 = 0,0,0


    def expansion(self, n):
        return (-1j)**n/(self.k*1j) * (2*n+1) / (n*(n+1));

    def BSC(self, n, m, mode='TE'):
        """Return the beam shape coefficients
         (:math:`g^{l}_{n, TE}`, :math:`g^{l}_{n, TM}`) for a plane wave.
         (Eq: VI.77 of G&G)
        Parameters
        ----------
        n : class:`int`
            Order of the expansion.
        m : class:`int`
            Description of parameter `m`.
        mode : class:`str`
            Mode of the plane wave, either 'TE' or 'TM'.
        Returns
        -------
        class:`float`
            Beam shape coefficient of order n and degree m.
        """
        if m not in [-1,1]: return 0

        if mode == 'TM': return 1 / 2

        if mode == 'TE': return  -m * 1j / 2

    def EField(self, x, y, z):
        """ Value of the electric Field
        :math:`E = E_0 \\exp{\big( i k (z-z_0) \big) (p_x \\vec{e_x}, p_y \\vec{e_y})}.

        """
        z0 = 0
        Pol = self.Polarization.Radian
        FPhase = exp(1j * self.k * (z - z0))
        return self.E0 * FPhase * cos(Pol), \
               self.E0 * FPhase * sin(Pol), \
               0

    def HField(self, x, y, z):
        """ Value of the electric Field
        :math:`H = H_0 \\exp{\big( i k (z-z_0) \big) (p_x \\vec{e_x}, p_y \\vec{e_y})}.

        """
        z0 = 0
        Pol = self.Polarization.Radian
        FPhase = exp(1j * self.k * (z - z0))
        return self.H0 * FPhase * cos(Pol), \
               self.H0 * FPhase * sin(Pol), \
               0




class GaussianBeam(BaseSource):
    def __init__(self,
                 Wavelength:   float,
                 NA:           float = 0.1,
                 Polarization: float = 0,
                 E0:           float = 1,
                 p:            int   = 20,
                 offset:       list = None):

        self.Wavelength = Wavelength
        self.k = 2 * pi / Wavelength
        self.Polarization = _Polarization(Polarization)
        self.E0 = E0
        self.NA = NA
        self.w0 = 2 * self.Wavelength / (pi * NA)
        self.s = 1/(self.k*self.w0)

        if offset is None:
            self.offset = array( [0,0,0] ) + 1e-16
        else:
            self.offset = array( offset )

        self.Offset = self.offset * self.k


        self.R0 = sqrt(self.Offset[0]**2 + self.Offset[1]**2)
        self.xi = arccos(self.Offset[0]/self.R0)

        termP = sqrt( 2.3 * p * ( self.s**(-2) + 4 * self.s**2 * self.Offset[2]**2 ) )

        NBSC = (max(1, self.R0 - 1/2 - termP), self.R0 - 1/2 + termP)


    def BSC(self, n, m, mode='TE'):
        """Return the beam shape coefficients
         (:math:`g^{l}_{n, TE}`, :math:`g^{l}_{n, TM}`) for a focused Gaussian
         beam using the quadrature method (ref[2]:Eq:17).

        Parameters
        ----------
        n : class:`int`
            Order of the expansion.
        m : class:`int`
            Description of parameter `m`.
        mode : class:`str`
            Mode of the plane wave, either 'TE' or 'TM'.
        Returns
        -------
        class:`float`
            Beam shape coefficient of order n and degree m.
        """
        if m not in [-1,1]: return 0

        if mode == 'TM': return self.Anm(n, m)

        if mode == 'TE': return self.Bnm(n, m)


    def expansion(self, n):
        return (-1j)**n/(self.k*1j) * (2*n+1) / (n*(n+1));


    def Phi0(self, x, y, z):
        """Method returns
        :math:`\\Phi_0 = -i Q \\exp \\Big[ i Q \\frac{(x-x_0)^2 + (y-y_0)^2}{w_0^2} \\Big]`

        """

        Q = self.Q(x,y,z)
        term = (x-self.offset[0])**2 + (y-self.offset[1])**2
        term *= 1j * Q / self.w0**2
        term = -1j * Q * exp(term)
        return term
        #return -1j * Q * exp(1j * Q * ( (x-self.offset[0])**2 + (y-self.offset[1])**2 )/self.w0**2 )


    def Phi0_Spherical(self, r, phi, theta):
        """Method returns
        :math:`\\Phi_0 = -i Q \\exp \\Big[ \\frac{i Q}{w_0^2} \\Big(  (r \\sin (\\theta) \\cos (\\phi) - x_0^2 ) + (r \\sin (\\theta) \\sin (\\phi) - y_0^2 ) \\Big) \\Big]`

        """
        Q = self.Q_Spherical(r, theta)
        term0 = -1j * Q
        term1 = 1j * Q / self.w0**2
        term2 = ( r * sin(theta) * cos(phi) - self.offset[0] )**2
        term3 = ( r * sin(theta) * sin(phi) - self.y0 )**2

        return term0 * exp(term1 * ( term2 + term3 ) )


    def Q(self,z):
        """Method returns
        :math:`Q_{cartesian} =\\frac{1}{2 ( z - z_0 )/k w_0^2 - i}`

        """
        return 1/( 2 * (z-self.offset[0])/(self.k * self.w0**2) - 1j )

    def Q_Spherical(self, r, theta):
        """Method returns
        :math:`Q_{spherical} = \\frac{1}{2 ( r \\cos (\\theta) - z_0 )/k w_0^2 - i}`

        """
        return 1/( 2 * ( r * cos(theta) -self.offset[0])/(self.k * self.w0**2 ) - 1j )

    def BSC(self, n, m, mode='TE'):
        """Return the beam shape coefficients
         (:math:`g^{l}_{n, TE}`, :math:`g^{l}_{n, TM}`) for a plane wave.
         (Eq: VI.77 of G&G)
        Parameters
        ----------
        n : class:`int`
            Order of the expansion.
        m : class:`int`
            Description of parameter `m`.
        mode : class:`str`
            Mode of the plane wave, either 'TE' or 'TM'.
        Returns
        -------
        class:`float`
            Expansion coefficient.
        """
        if m not in [-1,1]: return 0

        if mode == 'TM': return 1 / 2

        if mode == 'TE': return  -m * 1j / 2



    def Bnm(self, n, m, Num=200):
        """
        Method returns:

        :math:`B_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}\
         (\\beta) - \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_4 \\hat{I}_{m}\
          (\\beta) d \\theta`

        """
        Phi = linspace(0,2*pi,Num); Theta = linspace(0,pi,Num)

        rhon = (n + 0.5); r = rhon/self.k;

        sub = []
        for theta in Theta:
            Q = self.Q_Spherical(r, theta)
            beta = -2 * 1j * Q * self.s**2 * self.R0 * rhon *sin(theta)
            param = (n, m, rhon, Q, theta)
            term0 =  self.ImHat(m+1, beta) - self.ImHat(m-1, beta)
            term0 *= self.I_2(*param)
            term0 -= self.I_4(*param) * self.ImHat(m, beta)
            term0 *= self.I_1(*param)

            sub.append(-term0)

        sub = array(sub, copy=False).squeeze()

        return integrate( y=sub, x=Theta)* self.I_0(*param)


    def Anm(self, n, m, Num=200):
        """
        Method returns:

        :math:`A_{mn} = \\gamma_0 \\int_0^{\\pi} \\gamma_1 \\big[ \\hat{I}_{m+1}\
         (\\beta) + \\hat{I}_{m-1} (\\beta) \\big] -\\gamma_3 \\hat{I}_{m}\
          (\\beta) d \\theta`

        """
        Phi = linspace(0,2*pi,Num); Theta = linspace(0,pi,Num)

        rhon = (n + 0.5); r = rhon/self.k;

        sub = []
        for theta in Theta:
            Q = self.Q_Spherical(r, theta)
            beta = -2 * 1j * Q * self.s**2 * self.R0 * rhon *sin(theta)
            param = (n, m, rhon, Q, theta)
            term0 =  self.ImHat(m+1, beta) + self.ImHat(m-1, beta)
            term0 *= self.I_2(*param)
            term0 -= self.I_3(*param) * self.ImHat(m, beta)
            term0 *= self.I_1(*param)

            sub.append(-term0)

        sub = array(sub, copy=False).squeeze()

        return integrate( y=sub, x=Theta)* self.I_0(*param)


    def ImHat(self, m, beta):
        """
        Method returns:

        :math:`e^{-\\beta - i m \\xi} I_m(\\beta)`

        """
        return exp(-beta - 1j * m * self.xi ) * self.Im(m, beta)

    def Im(self, m, beta):
        """
        Method modified Bessel function :math:`I_m(\beta)` from ref[2]:Eq:16

        :math:`I_m(\\beta) = \\frac{1}{2\\pi} \int_0^{2\\pi}
        e^{z \\cos(\\phi) - i m \\phi} d\\phi`

        """
        x = linspace(0, 2*pi, 300)
        y = exp(beta * cos(x) - 1j*m*x)
        return 1/(2*pi) * integrate(y=y, x=x, axis=0)

    def I_0(self, n, m, rhon, Q, theta):
        """
        Method returns :math:`\gamma_0` from ref[2]:Eq:18

        :math:`\\frac{(-i)^n \\rho_n^2}{(2n+1) \\psi_n(\\rho_n)}e^{-i Z_0}`

        """
        term0 = (-1j)**n * (rhon)**2
        term1 = (2*n+1) * Xi(n, rhon)
        term2 = exp(-1j * self.Offset[2])

        return term0 / term1 * term2

    def I_1(self, n, m, rhon, Q, theta):
        """
        Method returns :math:`\gamma_1` from ref[2]:Eq:19

        :math:`Q \\exp \\big[ i Q s^2 \\Big( R_0 - \\rho_n \\sin (\\theta )  \\big)^2 \
        + i \\rho_n \\cos( \\theta ) \\Big] \\hat{P}_n^{|m|} (\\cos ( \\theta ) ) \
          \\sin(\\theta)`

        """
        factor = ( 2 * n + 1 ) * r8_factorial(n - m)
        factor /= ( 2 * r8_factorial( n + m ) )

        term1 = 1j * Q * self.s**2
        term1 *= ( self.R0 - rhon * sin(theta) )**2
        term3 = 1j * rhon * cos(theta)
        term4 = sqrt(factor) * Pnm_(n, abs(m), cos([theta]))[-1] * sin(theta)

        return Q * exp(term1 + term3 ) * term4

    def I_2(self, n, m, rhon, Q, theta):
        """
        Method returns :math:`\gamma_2` from ref[2]:Eq:19

        :math:`\\big( 2 Q s^2 \\rho_n \\cos(\\theta) -1 \\big) \\sin(\\theta)`

        """
        return ( 2 * Q * self.s**2 * rhon * cos(theta) - 1 ) * sin(theta)

    def I_3(self, n, m, rhon, Q, theta):
        """
        Method returns :math:`\gamma_3` from ref[2]:Eq:19

        :math:`4 Q s^2 X_0 \\cos(\\theta)`

        """
        return 4 * Q * self.s**2 * self.Offset[0] * cos(theta)

    def I_4(self, n, m, rhon, Q, theta):
        """
        Method returns :math:`\gamma_4` from ref[2]:Eq:19

        :math:`4 Q s^2 Y_0 \\cos(\\theta)`

        """
        return 4 * Q * self.s**2 * self.Offset[1] * cos(theta)


    def EField(self, x, y, z):
        """ Value of the electric Field (First Order)
        ref[2]:Eq:6
        :math:`E = E_0 \\exp{\big( i k (z-z_0) \big) (p_x \\vec{e_x}, p_y \\vec{e_y})}`.

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
        """ Value of the electric Field (First Order)
        ref[2]:Eq:6
        :math:`H = H_0 \\exp{\big( i k (z-z_0) \big) (p_x \\vec{e_x}, p_y \\vec{e_y})}`.

        """

        s = 1 / (self.k * self.w0)
        _Q = self.Q(x,y,z)
        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * (z - self.offset[2]))

        return 0, \
               self.H0 * self.Phi0(x,y,z) * FPhase * 1, \
               self.H0 * self.Phi0(x,y,z) * FPhase * -2*s*_Q*(y-self.y0)/self.w0


    def Er(self, r, phi, theta):
        """ Value of the electric Field (First Order)
        ref[2]:Eq:10
        :math:`E_r = E_0 \\Phi_0 \\Big[ \\sin (\\theta) \\cos (\\phi) - \\frac{2Q}{kw_0^2} \\cos (\\theta) \\big( r \\sin (\\theta) \\cos (\\phi) -x_0 \\big)  \\Big] \\exp^{ik (r \\cos (\\theta) -z_0 )}`.

        """

        w0 = 1e-5

        s = 1 / (self.k * w0)

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
        """ Value of the magnetic Field (First Order)
        ref[2]:Eq:10
        :math:`H_r = E_0 \\Phi_0 \\Big[ \\sin (\\theta) \\cos (\\phi) - \\frac{2Q}{kw_0^2} \\cos (\\theta) \\big( r \\sin (\\theta) \\cos (\\phi) -y_0 \\big)  \\Big] \\exp^{ik (r \\cos (\\theta) -z_0 )}`.

        """


        s = 1 / (self.k * self.w0)

        _Q = self.Q_Spherical(r, theta)

        Pol = self.Polarization.Radian

        FPhase = exp(1j * self.k * ( r * cos(theta) - self.offset[2] ) )

        term0 = self.H0 * self.Phi0_Spherical(r, phi, theta)

        term1 = sin(theta) * cos(phi)

        term2 =  2 * _Q / ( self.k * self.w0**2 ) * cos(theta)

        term3 = (r * sin(theta) * cos(phi) - self.y0)

        Hr = term0 * (term1 - term2 * term3) * FPhase

        return Hr
