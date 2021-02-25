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





def NPnm( n, m, x ):
    #ref : https://people.sc.fsu.edu/~jburkardt/py_src/polpak/legendre_associated_normalized.py


  if ( m < 0 ):
    print ( '' )
    print ( 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!' )
    print ( '  Input value of M is %d' % ( m ) )
    print ( '  but M must be nonnegative.' )

  if ( n < m ):
    print ( '' )
    print ( 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!' )
    print ( '  Input value of M = %d' % ( m ) )
    print ( '  Input value of N = %d' % ( n ) )
    print ( '  but M must be less than or equal to N.' )

  if ( x < -1.0 ):
    print ( '' )
    print ( 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!' )
    print ( '  Input value of X = %f' % ( x ) )
    print ( '  but X must be no less than -1.' )

  if ( 1.0 < x ):
    print ( '' )
    print ( 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!' )
    print ( '  Input value of X = %f' % ( x ) )
    print ( '  but X must be no more than 1.' )

  cx = np.zeros ( n + 1 )

  cx[m] = 1.0
  somx2 = np.sqrt ( 1.0 - x * x )

  fact = 1.0
  for i in range ( 0, m ):
    cx[m] = - cx[m] * fact * somx2
    fact = fact + 2.0

  if ( m != n ):

    cx[m+1] = x * float ( 2 * m + 1 ) * cx[m]

    for i in range ( m + 2, n + 1 ):
      cx[i] = ( float ( 2 * i     - 1 ) * x * cx[i-1] \
              + float (   - i - m + 1 ) *     cx[i-2] ) \
              / float (     i - m     )

  for mm in range ( m, n + 1 ):

    factor = np.sqrt ( ( ( 2 * mm + 1 ) * r8_factorial ( mm - m ) ) \
      / (  r8_factorial ( mm + m ) ) )

    cx[mm] = cx[mm] * factor

  return cx

def r8_factorial ( n ):


  from sys import exit

  if ( n < 0 ):
    print ( '' )
    print ( 'R8_FACTORIAL - Fatal error!' )
    print ( '  N < 0.' )
    exit ( 'R8_FACTORIAL - Fatal error!' )

  value = 1.0

  for i in range ( 2, n + 1 ):
    value = value * i

  return value
#-
