import numpy as np
from scipy import special
from numpy import cos, sin
from scipy.special import spherical_yn, spherical_jn, hankel1

def riccati1(n,x):
    jn = spherical_jn(n, x)
    jnp = spherical_jn(n, x, derivative=True)

    return np.array([x*jn, jn + x*jnp])

def riccati1_pp(n, x):
    return 1 / x * (n + n**2 - x**2) * spherical_jn(n, x)

def _riccati_bessel_j(n, x):
    return special.riccati_jn(n, x)


def riccati2(n,x):
    yn = spherical_yn(n, x)
    ynp = spherical_yn(n, x, derivative=True)

    return np.array([x*yn, yn + x*ynp])


def riccati3(n,x):
    return riccati_1_single(n,x) + 1j*riccati_2_single(n,x)

def spherical_hn(n, z, derivative=False):
    return spherical_jn(n,z,derivative) + 1j*spherical_yn(n,z,derivative)


def Pin(n, theta):
    lpn = special.legendre(n)
    lpn_p = lpn.deriv()

    return -1*lpn_p(cos(theta))


def Taun(n, theta):
    lpn = special.legendre(n)
    lpn_p = lpn.deriv()
    lpn_p2 = lpn_p.deriv()

    return -1*cos(theta)*lpn_p(cos(theta)) + sin(theta)**2*lpn_p2(cos(theta))



def M1o1n(n, k, r, theta, phi):
    theta_comp =      cos(phi) * Pin(n, theta)  * spherical_jn(n, k * r)
    phi_comp   = -1 * sin(phi) * Taun(n, theta) * spherical_jn(n, k * r)
    r_comp     = np.zeros(shape = theta.shape, dtype=np.complex)

    return np.array([r_comp, theta_comp, phi_comp])


def M1e1n(n, k, r, theta, phi):
    theta_comp = -1*sin(phi) * Pin(n, theta)  * spherical_jn(n, k * r)
    phi_comp =   -1*cos(phi) * Taun(n, theta) * spherical_jn(n, k * r)
    r_comp = np.zeros(shape = theta.shape, dtype=np.complex)

    return np.array([r_comp, theta_comp, phi_comp])


def N1o1n(n, k, r, theta, phi):
    p = k * r
    theta_comp = sin(phi) * Taun(n, theta) * (spherical_jn(n, p) + p * spherical_jn(n, p, derivative=True) ) / p
    phi_comp   = cos(phi) * Pin(n, theta)  * (spherical_jn(n, p) + p * spherical_jn(n, p, derivative=True) ) / p
    r_comp     = sin(phi) * n * (n + 1) * sin(theta) * Pin(n, theta) * spherical_jn(n, p) / p

    return np.array([r_comp, theta_comp, phi_comp])



def N1e1n(n, k, r, theta, phi):
    p = k*r
    theta_comp =    cos(phi) * Taun(n, theta) * ( spherical_jn(n, p) + p * spherical_jn(n, p, derivative=True) ) / p
    phi_comp   = -1*sin(phi) * Pin(n, theta)  * ( spherical_jn(n, p) + p * spherical_jn(n, p, derivative=True) ) / p
    r_comp     =    cos(phi) * n * (n+1)      * sin(theta) * Pin(n, theta) * spherical_jn(n, p) / p

    return np.array([r_comp, theta_comp, phi_comp])




def M3o1n(n, k, r, theta, phi):
    theta_comp =      cos(phi) * Pin(n, theta)  * spherical_hn(n, k * r)
    phi_comp   = -1 * sin(phi) * Taun(n, theta) * spherical_hn(n, k * r)
    r_comp     = np.zeros(shape = theta.shape, dtype=np.complex)

    return np.array([r_comp, theta_comp, phi_comp])


def M3e1n(n, k, r, theta, phi):
    theta_comp = -1*sin(phi) * Pin(n, theta)  * spherical_hn(n, k * r)
    phi_comp =   -1*cos(phi) * Taun(n, theta) * spherical_hn(n, k * r)
    r_comp = np.zeros(shape = theta.shape, dtype=np.complex)

    return np.array([r_comp, theta_comp, phi_comp])


def N3o1n(n, k, r, theta, phi):
    p = k * r
    theta_comp = sin(phi) * Taun(n, theta) * (spherical_hn(n, p) + p * spherical_hn(n, p, derivative=True) ) / p
    phi_comp   = cos(phi) * Pin(n, theta)  * (spherical_hn(n, p) + p * spherical_hn(n, p, derivative=True) ) / p
    r_comp     = sin(phi) * n * (n + 1) * sin(theta) * Pin(n, theta) * spherical_hn(n, p) / p

    return np.array([r_comp, theta_comp, phi_comp])



def N3e1n(n, k, r, theta, phi):
    p = k*r
    theta_comp =    cos(phi) * Taun(n, theta) * ( spherical_hn(n, p) + p * spherical_hn(n, p, derivative=True) ) / p
    phi_comp   = -1*sin(phi) * Pin(n, theta)  * ( spherical_hn(n, p) + p * spherical_hn(n, p, derivative=True) ) / p
    r_comp     =    cos(phi) * n * (n+1)      * sin(theta) * Pin(n, theta) * spherical_hn(n, p) / p

    return np.array([r_comp, theta_comp, phi_comp])



def mie_coefficient_a(order,
                      diameter=10E-6,
                      permeability=1.,
                      sp_permeability=1.,
                      wavelength=1e-6,
                      reffractive=1.4):
    """ Obtains the value of Mie coefficient :math:`a_n`.
    Due to the magnitude of the denominator for sufficiently high orders,
    we treat the case where this denominator is too high or zero - which
    happens when its value is too high to be represented in a computer.
    """
    alpha = np.pi * diameter / wavelength
    beta = reffractive * alpha
    denominator = (sp_permeability * riccati_bessel_y(order, alpha)
                   * d_riccati_bessel_j(order, beta)
                   - permeability * reffractive
                   * d_riccati_bessel_y(order, alpha)
                   * riccati_bessel_j(order, beta))
    if denominator > 1E20 or not denominator:
        return 0

    return (sp_permeability * riccati_bessel_j(order, alpha)
            * d_riccati_bessel_j(order, beta)
            - permeability * reffractive
            * d_riccati_bessel_j(order, alpha)
            * riccati_bessel_j(order, beta)) / denominator

def expansion(order, k):
    return (-1j)**order/(k*1j) * (2*order+1) / (order*(order+1));



def riccati_bessel_y(degree, argument):
    """ Riccati-Bessel function of second kind """
    return special.riccati_jn(degree, float(argument))[0][degree] \
           - 1j * special.riccati_yn(degree, float(argument))[0][degree]


def d_riccati_bessel_y(degree, argument):
    """ Derivative of Riccati-Bessel function of second kind """
    return special.riccati_jn(degree, float(argument))[1][degree] \
           - 1j * special.riccati_yn(degree, float(argument))[1][degree]


def riccati_bessel_j(degree, argument):
    """ Riccati-Bessel function of first kind """
    return special.riccati_jn(degree, float(argument))[0][degree]

def d_riccati_bessel_j(degree, argument):
    """ Derivative of Riccati-Bessel function of first kind """
    return special.riccati_jn(degree, float(argument))[1][degree]

from ai import cs
from mayavi import mlab


def legendre_p(degree, order, argument):
    if degree < np.abs(order):
        return 0
    if order < 0:
        return pow(-1, -order) / fac_plus_minus(degree, -order) * legendre_p(degree, -order, argument)
    return special.lpmv(order, degree, argument)

def _riccati_bessel_j(degree, argument):
    return special.riccati_jn(degree, argument)

def _riccati_bessel_y(degree, argument):
    return (np.array(special.riccati_jn(degree, argument))
            - 1j * np.array(special.riccati_yn(degree, argument)))

w = 1e-6; k = 2*np.pi/w;
PHI, THETA = np.mgrid[-np.pi/2:np.pi/2:complex(0,100),-np.pi:np.pi:complex(0,100)]
R = PHI*0+1e-4


def get_max_it(a):
    return 4

def d2_riccati_bessel_y(degree, argument):

    if argument == 0:
        argument = 1E-16
    return 1 / (argument) \
            * (degree + pow(degree, 2) - pow(argument, 2)) \
            * (special.spherical_jn(degree, float(argument)) \
               - 1j * special.spherical_yn(degree, float(argument)))







N = 5
coef = mie_coefficient_a(N)

for n in range(1, N):
    for l in [-1,1]:

        increment = coef(n) *








"""


result = np.array(result).reshape(THETA.shape)

r, t, p = N3e1n(n=1, k=k,r=R, theta=THETA, phi = PHI)

X, Y, Z = cs.sp2cart(result, PHI, THETA)


mlab.mesh(X.real, Y.real, Z.real)
mlab.show()
"""

#-
