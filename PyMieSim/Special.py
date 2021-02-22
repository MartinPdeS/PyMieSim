import numpy as np
from scipy import special
from numpy import cos, sin, sqrt
from scipy.special import spherical_yn as yn, spherical_jn as jn, hankel1 as h
from scipy.special import legendre as Pn, factorial as fac, lpmn, gammasgn, loggamma



def hn(n, z, derivative=False):
    return jn(n,z,derivative) + 1j*yn(n,z,derivative)


def mnFactorial(n, m):
    """ Calculates the expression below avoiding overflows.

    .. math::
        \\frac{(n + m)!}{(n - m)!}
    """
    return gammasgn(n) * gammasgn(m) * exp(loggamma(n) - loggamma(m))


@np.vectorize
def Pnm(n, m, x):
    """Eq:II.77 """
    return lpmn(m, n, x)[0][-1, -1]


@np.vectorize
def Pnm_p(n, m, x):
    """Eq:II.77 """
    return lpmn(m, n, x)[1][-1, -1]


def Taun(n, x):
    return -1*Pn(n).deriv()(cos(x))


def Pin(n, x):
    return -1*cos(x)*Pn(n).deriv(1)(cos(x)) + sin(x)**2*Pn(n).deriv(2)(cos(x))


def Taunm(n, m, x):
    """Eq: III.51 """
    index = np.isclose(x,1,1e-6)
    if len(index) > 0: x[index] = 1-1e-6
    index = np.isclose(x,-1,1e-6)
    if len(index) > 0: x[index] = -1+1e-6

    return sqrt(1-x**2) * Pnm_p(n, m, x)


def Pinm(n, k, x):
    """Eq: III.52 """
    index = np.isclose(x,1,1e-6)
    if len(index) > 0: x[index] = 1-1e-6
    index = np.isclose(x,-1,1e-6)
    if len(index) > 0: x[index] = -1+1e-6

    return -Pnm(n, k, x) / (sqrt(1-x**2))


def _Psi(type, n, x):
    """Eq:II.83-86 """
    if type == 0: return x*jn(n, x)
    if type == 1: return jn(n, x)
    if type == 2: return yn(n, x)
    if type == 3: return _Psi(1, n, x) + 1j * _Psi(2, n, x)
    if type == 4: return _Psi(1, n, x) - 1j * _Psi(2, n, x)


def Psi(n, x):
    return x * _Psi(1, n, x)

def Psi_p(n, x):
    return x * _Psi_p(1, n, x) + _Psi(1, n, x)


def _Psi_p(type, n, x):
    """Eq:II.83-86 """
    if type == 0: return x*jn(n, x, derivative=True) + jn(n, x, derivative=False)
    if type == 1: return jn(n, x, derivative=True)
    if type == 2: return yn(n, x, derivative=True)
    if type == 3: return jn(n,x,derivative=True) + 1j*yn(n,x,derivative=True)
    if type == 4: return jn(n,x,derivative=True) - 1j*yn(n,x,derivative=True)


def Xi(n, x):
    """Eq:II.87 """
    return x * _Psi(4, n, x)

def Xi_p(n, x):
    return x * _Psi_p(4, n, x) + _Psi(4, n, x)








class GeneralFunc():
    def __init__(self, function, *args):
        self.functions = [function]
        self.Operation = ['.']
        self.args = [args]

    def eval(self):

        res = complex(0)
        for n, func in enumerate(self.functions):
            Op = self.Operation[n]; arg = self.args[n]
            if arg == '':
                if Op ==  '.': res += func; print('Identity')
                if Op ==  '+': res += func; print('Addition')
                if Op ==  '-': res -= func; print('Subtraction')
                if Op ==  '*': res *= func; print('Multiplication')
                if Op ==  '/': res /= func; print('Division')
                if Op == '**': res = res ** func; print('Division')
            else:
                if Op ==  '.': res += func(*args); print('Identity')
                if Op ==  '+': res += func(*args); print('Addition')
                if Op ==  '-': res -= func(*args); print('Subtraction')
                if Op ==  '*': res *= func(*args); print('Multiplication')
                if Op ==  '/': res /= func(*args); print('Division')
                if Op == '**': res = res ** func(*args); print('Division')

        return res


    def __add__(self, func):
        self.Operation.append('+')
        self.functions.append(func.functions[0])
        self.args.append(func.args[0])


    def __sub__(self, func):
        self.Operation.append('-')
        self.functions.append(func.functions[0])
        self.args.append(func.args[0])

    def __mul__(self, func):
        self.Operation.append('*')
        self.functions.append(func.functions[0])
        self.args.append(func.args[0])

    def __div__(self, func):
        self.Operation.append('/')
        self.functions.append(func.functions[0])
        self.args.append(func.args[0])

    def __pow__(self, func):
        self.Operation.append('**')

        if not isinstance(func, list):
            self.functions.append(func)
            self.args.append('')
        else:
            self.functions.append(func.functions[0])
            self.args.append(func.args[0])

#-
