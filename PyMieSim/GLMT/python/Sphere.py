# -*- coding: utf-8 -*-
# http://pymiescatt.readthedocs.io/en/latest/forward.html

from numpy import zeros, cos, outer, exp, pi, ndarray

from PyMieSim.Tools.Special import Taunm, Pinm
from PyMieSim.BaseClasses import BaseScatterer, BaseSource

def GetMaxOrder(Scat):
    """Function compute max Order (n) for :math:`S_1 \, \\& \, S_2` computing
    MaxOrder = :math:`\\big( 2+x + 4 x^{1/3} \\big)`;

    x begin the size parameter of the scatterer.

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.

    Returns
    -------
    :class:`int`
        The maximum order n to reach.

    """
    return int(2 + Scat.SizeParam + 4*Scat.SizeParam**(1/3))

def S1(Scat,
       Source,
       phi,
       theta):
    """Function compute :math:`S_2` component of scattered FarField defined
    as Eq:III.110 of B&B.
    :math:`S_2 = \sum_{n=1}^\\infty \sum_{m=-n}^{n} \\frac{2n+1}{n(n+1)}
    \\Big[ m a_n g_{n,TM}^m  \\pi_n^{|m|} \\big(\cos (\\theta) \\big) +
    i b_n g_{n,TE}^m \\tau_n^{|m|} \\big(\\cos (\\theta) \\big) \\Big] \\exp{\\big(i m \phi\\big)}`

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.
    Source : :class:`BaseSource`
        Description of parameter `Source`.
    phi : :class:`numpy.ndarray`
        Vector representing :math:`\\phi` space :math:`[-\\frac{\pi}{2}: \\frac{\pi}{2}]`.
    theta : :class:`numpy.ndarray`
        Vector representing :math:`\\theta` space :math:`[-\pi: \pi]`.

    Returns
    -------
    :class:`numpy.ndarray`
        2D ndarray representing :math:`S_1\\big(\\phi, \\theta \\big)`.

    """
    result = zeros([phi.size, theta.size], dtype=complex)
    MaxOrder = GetMaxOrder(Scat)
    #phi -= pi/2  # conversion from spherical coordinates to "LMT" spherical
    an = Scat.an(MaxOrder); bn = Scat.bn(MaxOrder)
    for n in range(1, MaxOrder):
        prefactor = (2*n+1)/(n*(n+1))

        for m in [-1,1]:

            Lterm =  m * an[n] * Source.BSC(n, m, 'TM') * Pinm(n,abs(m), cos(phi-pi/2))

            Rterm = 1j * bn[n] * Source.BSC(n, m, 'TE') * Taunm(n,abs(m), cos(phi-pi/2))

            result += outer(prefactor*( Rterm  + Lterm ), exp(1j*m*theta))

    return result


def S2(Scat,
       Source,
       phi,
       theta):
    """Function compute :math:`S_2` component of scattered FarField defined
    as Eq:III.111 of B&B.
    :math:`S_2 = \sum_{n=1}^\\infty \sum_{m=-n}^{n} \\frac{2n+1}{n(n+1)}
    \\Big[ a_n g_{n,TM}^m  \\tau_n^{|m|} \\big(\cos (\\theta) \\big) +
    m i b_n g_{n,TE}^m \\pi_n^{|m|} \\big(\\cos (\\theta) \\big) \\Big] \\exp{\\big(i m \phi\\big)}`

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.
    Source : :class:`BaseSource`
        Description of parameter `Source`.
    phi : :class:`numpy.ndarray`
        Vector representing :math:`\\phi` space :math:`[-\\frac{\pi}{2}: \\frac{\pi}{2}]`.
    theta : :class:`numpy.ndarray`
        Vector representing :math:`\\theta` space :math:`[-\pi: \pi]`.

    Returns
    -------
    :class:`numpy.ndarray`
        2D ndarray representing :math:`S_1\\big(\\phi, \\theta \\big)`.

    """
    result = zeros([phi.size, theta.size], dtype=complex)
    MaxOrder = GetMaxOrder(Scat)
    #phi -= pi/2  # conversion from spherical coordinates to "LMT" spherical

    an = Scat.an(MaxOrder); bn = Scat.bn(MaxOrder)
    for n in range(1, MaxOrder):
        prefactor = (2*n+1)/(n*(n+1))

        for m in [-1,1]:

            Lterm =  an[n] * Source.BSC(n, m, 'TM') * Taunm(n,abs(m), cos(phi-pi/2))

            Rterm =  1j * m * bn[n] * Source.BSC(n, m, 'TE') * Pinm(n,abs(m), cos(phi-pi/2))

            result += outer(prefactor*( Rterm  + Lterm ), exp(1j*m*theta))

    return result


def SPF(Scat,
        Source,
        phi:     ndarray,
        theta:   ndarray,
        r:       float = 1):

    """Function compute the scattering phase function (SPF) defined as
    :math:`SPF = \\sqrt{ |E_\\theta |^2 + |E_\\phi |^2 }`

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.
    Source : :class:`BaseSource`
        Description of parameter `Source`.
    phi : ndarray
        Vector representing :math:`\\phi` space :math:`[-\\frac{\pi}{2}: \\frac{\pi}{2}]`.
    theta : ndarray
        Vector representing :math:`\\theta` space :math:`[-\pi: \pi]`.
    r : float
        Distance of the FarField


    Returns
    -------
    :class:`numpy.ndarray`
        2D ndarray representing SPF.

    """

    EPhi = FieldPhi(Scat, Source, phi, theta, r)

    ETheta = FieldTheta(Scat, Source, phi, theta, r)

    return abs(EPhi)**2 + abs(ETheta)**2


def FieldTheta(Scat,
               Source,
               phi:     ndarray,
               theta:   ndarray,
               r:       float =1):
    """Function compute :math:`\vec{\\theta}` component of scattered FarField defined
    as Eq:III.108 of B&B.
    :math:`E_{\\theta} = \\frac{i E_0}{kr} \exp{(-ikr)} S_2`

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.
    Source : :class:`BaseSource`
        Description of parameter `Source`.
    phi : :class:`numpy.ndarray`
        Vector representing :math:`\\phi` space :math:`[-\\frac{\pi}{2}: \\frac{\pi}{2}]`.
    theta : :class:`numpy.ndarray`
        Vector representing :math:`\\theta` space :math:`[-\pi: \pi]`.
    r : :class:`float`
        Distance of the FarField

    Returns
    -------
    :class:`numpy.ndarray`
        2D ndarray representing the scattered field.

    """


    ETheta = 1j * Source.E0 / (Source.k*r) \
           * exp(-1j*Source.k*r)    \
           * S1(Scat, Source, phi, theta)

    return ETheta



def FieldPhi(Scat,
             Source,
             phi:     ndarray,
             theta:   ndarray,
             r:       float = 1):
    """Function compute :math:`\\vec{\\phi}` component of scattered FarField defined
    as Eq:III.109 of B&B.
    :math:`E_{\\theta} = -\\frac{E_0}{kr} \\exp{(-ikr)} S_1`

    Parameters
    ----------
    Scat : :class:`BaseScatterer`
        Description of parameter `Scat`.
    Source : :class:`BaseSource`
        Description of parameter `Source`.
    phi : :class:`numpy.ndarray`
        Vector representing :math:`\\phi` space :math:`[-\\frac{\pi}{2}: \\frac{\pi}{2}]`.
    theta : :class:`numpy.ndarray`
        Vector representing :math:`\\theta` space :math:`[-\pi: \pi]`.
    r : :class:`float`
        Distance of the FarField

    Returns
    -------
    :class:`numpy.ndarray`
        2D ndarray representing the scattered field.

    """

    EPhi = - Source.E0 / (Source.k*r) \
          * exp(-1j*Source.k*r)    \
          * S2(Scat, Source, phi, theta)

    return EPhi
