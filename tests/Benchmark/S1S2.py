import matplotlib.pyplot as plt
import numpy as np
import timeit

from PyMieSim.LMT.python.Sphere import S1S2 as PyS1S2
from PyMieSim.LMT.Scatterer import SPHERE


def Speed(setup):

    BenchPython = """PyS1S2(*args);"""

    BenchCpp = """scat.S1S2(Phi);"""

    print('\nPYTHON BENCHMARK')

    Bench = timeit.timeit(setup = setup,stmt = BenchPython, number = 1)

    print(Bench)

    print('C++ BENCHMARK\n' + '='*50)

    Bench = timeit.timeit(setup = setup,stmt = BenchCpp, number = 1)

    print(Bench)




def Correctness():
    Phi = np.linspace(-np.pi,np.pi,100);

    args = (1.8, 3e-6, 1e-6, 1, Phi)

    PyRes = PyS1S2(*args)

    scatt = SPHERE(Index = 1.8, Diameter=3e-6,Wavelength=1e-6)


    CppS1, CppS2 = scatt.S1S2(Phi+np.pi/2)

    print('##########', CppS2.shape)

    fig = plt.figure(figsize=(10,5))

    ax0 = fig.add_subplot(211); ax1 = fig.add_subplot(212)

    ax0.plot(Phi, np.imag(PyRes[0]), 'C0', linewidth=2, label='Python S1');
    ax1.plot(Phi, np.imag(PyRes[1]), 'C0',linewidth=2, label='Python S2')

    ax1.plot(Phi, np.imag(CppS2), 'C1--',linewidth=2, label='C++ S1');
    ax0.plot(Phi, np.imag(CppS1), 'C1--',linewidth=2, label='C++ S2')

    ax0.grid(); ax1.grid();
    ax0.legend(); ax1.legend();

    plt.show()





setup = """
import numpy as np
from PyMieSim.LMT.python.Sphere import S1S2 as PyS1S2
from PyMieSim.LMT.Scatterer import SPHERE
Phi = np.linspace(-np.pi,np.pi,20000)
args = (1.8, 3e-6, 1e-6, 1, Phi)
scat = SPHERE(Index = 1.8, Diameter=3e-6,Wavelength=1e-6)
"""

if __name__ == '__main__':

    Speed(setup)
    Correctness()










    # -
