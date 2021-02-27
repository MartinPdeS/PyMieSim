import matplotlib.pyplot as plt
import numpy as np
import timeit

from PyMieSim.LMT.python.Sphere import S1S2 as PyS1S2
from PyMieSim.LMT.Sphere import S1S2 as CppS1S2


def Speed(setup):

    BenchPython = """PyS1S2(*args);"""

    BenchCpp = """CppS1S2(*args);"""

    print('\nPYTHON BENCHMARK')

    Bench = timeit.timeit(setup = setup,stmt = BenchPython, number = 1)

    print(Bench)

    print('C++ BENCHMARK\n' + '='*50)

    Bench = timeit.timeit(setup = setup,stmt = BenchCpp, number = 1)

    print(Bench)




def Correctness():
    Phi = np.linspace(-np.pi,np.pi,100);

    args = (1.8, 3e-6, 1e-6, 1, Phi)
    cppRes = CppS1S2(*args)
    PyRes = PyS1S2(*args)


    fig = plt.figure(figsize=(10,5))

    ax0 = fig.add_subplot(211); ax1 = fig.add_subplot(212)


    ax0.plot(Phi, np.abs(PyRes[0]), 'C0', linewidth=2, label='Python S1');
    ax1.plot(Phi, np.abs(PyRes[1]), 'C0',linewidth=2, label='Python S2')

    ax0.plot(Phi, np.abs(cppRes[0]), 'C1--',linewidth=2, label='C++ S1');
    ax1.plot(Phi, np.abs(cppRes[1]), 'C1--',linewidth=2, label='C++ S2')

    ax0.grid(); ax1.grid();
    ax0.legend(); ax1.legend();

    plt.show()





setup = """
import numpy as np
from PyMieSim.LMT.python.Sphere import S1S2 as PyS1S2
from PyMieSim.LMT.Sphere import S1S2 as CppS1S2
Phi = np.linspace(-np.pi,np.pi,10000)
args = (1.8, 3e-6, 1e-6, 1, Phi)
"""

if __name__ == '__main__':

    Speed(setup)
    Correctness()










    # -
