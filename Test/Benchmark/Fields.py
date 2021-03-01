import matplotlib.pyplot as plt
import numpy as np
import timeit

from PyMieSim.utils import PlotField
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.LMT.Sphere import Fields as PyBindFields




def Speed(setup):
    BenchPython = """PyField(*args0)"""

    BenchCpp    = """CppField(*args0)"""

    BenchPyBind = """PyBindFields(*args1)"""

    print('\nPYTHON BENCHMARK')

    Bench = timeit.timeit(setup = setup,stmt = BenchPython, number = 1)

    print(Bench)

    print('='*50 + '\C++ BENCHMARK')

    Bench = timeit.timeit(setup = setup,stmt = BenchPyBind, number = 1)

    print(Bench)



def Correctness():
    Phi = np.linspace(-np.pi/2,np.pi/2,100); Theta = np.linspace(-np.pi,np.pi,120)
    THETA, PHI = np.meshgrid(Theta, Phi)

    args0 = (1.4, 10e-6, 1e-6, 1, Phi.flatten(), Theta.flatten(), 0,1,1)
    args1 = (1.4, 10e-6, 1e-6, 1, PHI.flatten(), THETA.flatten(), 0,1,1, THETA.flatten().size)

    PyParallel, PyPerpendicular = PyField(*args0);

    CppParallel, CppPerpendicular = PyBindFields(*args1);

    CppParallel = np.array(CppParallel).reshape([Phi.size, Theta.size])
    CppPerpendicular = np.array(CppPerpendicular).reshape([Phi.size, Theta.size])

    PlotField(Theta, Phi, CppParallel, CppPerpendicular)

    PlotField(Theta, Phi, PyParallel, PyPerpendicular)

    plt.show()



setup = """
import numpy as np
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.LMT.Sphere import Fields as PyBindFields
Phi = np.linspace(-np.pi/2,np.pi/2,200); Theta = np.linspace(-np.pi,np.pi,520)
PHI, THETA = np.meshgrid(Theta, Phi)
args = (1.8, 3e-6, 1e-6, 1, Phi)
args0 = (1.4, 1e-6, 1e-6, 1, Phi.flatten(), Theta.flatten(), 0,1,1)
args1 = (1.4, 1e-6, 1e-6, 1, PHI.flatten(), THETA.flatten(), 0,1,1, Phi.flatten().size)

"""

if __name__=='__main__':
    Speed(setup)
    Correctness()













    # -
