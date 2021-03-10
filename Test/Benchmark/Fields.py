import matplotlib.pyplot as plt
import numpy as np
import timeit

from PyMieSim.utils import PlotField
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.GLMT.Sphere import FieldsStructured as GLMTCppFields
from PyMieSim.LMT.Sphere import FieldsStructured as LMTCppFields




def Speed(setup):
    BenchPython = """PyField(*args0)"""

    LMTBenchPyBind = """LMTCppFields(*args1)"""

    GLMTLMTBenchPyBind = """GLMTCppFields(*args2)"""

    Bench = timeit.timeit(setup = setup,stmt = BenchPython, number = 1)

    print('\nLMT PYTHON BENCHMARK: ', Bench)

    Bench = timeit.timeit(setup = setup,stmt = LMTBenchPyBind, number = 1)

    print('\n\n'+'='*50 + '\n\n LMT C++ BENCHMARK', Bench)

    Bench = timeit.timeit(setup = setup,stmt = GLMTLMTBenchPyBind, number = 1)

    print('\n\n'+'='*50 + '\n\n GLMT C++ BENCHMARK', Bench)


def Correctness():
    Phi = np.linspace(-np.pi/2,np.pi/2,100); Theta = np.linspace(-np.pi,np.pi,100)
    THETA, PHI = np.meshgrid(Theta, Phi)

    args0 = (1.4, 10e-6, 1e-6, 1, Phi, Theta, 0,1,1)
    args1 = (1.4, 10e-6, 1e-6, 1, Phi, Theta, 0,1,1)

    PyParallel, PyPerpendicular = PyField(*args0);

    CppParallel, CppPerpendicular = LMTCppFields(*args1);

    PlotField(Theta, Phi, CppParallel, CppPerpendicular)

    PlotField(Theta, Phi, PyParallel, PyPerpendicular)

    plt.show()



setup = """
import numpy as np
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.LMT.Sphere import FieldsStructured as LMTCppFields
from PyMieSim.GLMT.Sphere import FieldsStructured as GLMTCppFields
from PyMieSim.Source import PlaneWave

beam = PlaneWave(Wavelength=1e-6)
BSC = beam.GetBSC(MaxOrder=10)
Phi = np.linspace(-np.pi/2,np.pi/2,100); Theta = np.linspace(-np.pi,np.pi,100)
PHI, THETA = np.meshgrid(Theta, Phi)

args0 = (1.4, 1e-6, 1e-6, 1, Phi.flatten(), Theta.flatten(), 0,1,1)
args1 = (1.4, 1e-6, 1e-6, 1, Phi.flatten(), Theta.flatten(), 0,1,1)
args2 = (1.4, 1e-6, 1e-6, 1, Phi, Theta, 0,1,1, beam._BSC_, beam.MaxOrder)



"""

if __name__=='__main__':
    Speed(setup)
    Correctness()


"""

0.16491000000678468

0.130888280000363
"""











    # -
