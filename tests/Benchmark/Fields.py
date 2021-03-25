import matplotlib.pyplot as plt
import numpy as np
import timeit

from PyMieSim.Plots import PlotField
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.Source import PlaneWave
from PyMieSim.LMT.Scatterer import SPHERE as SPHERELMT
from PyMieSim.GLMT.Scatterer import SPHERE as SPHEREGLMT


def Speed(setup):
    BenchPython = """PyField(**argsPy)"""

    GLMTClass = """B.sFields(Phi=Phi, Theta=Theta, R=1. );"""

    LMTClass = """A.sFields(Phi=Phi, Theta=Theta, R=1. ); """

    Bench = timeit.timeit(setup = setup,stmt = BenchPython, number = 1)
    print('\nLMT PYTHON BENCHMARK: ', Bench)

    Bench = timeit.timeit(setup = setup,stmt = LMTClass, number = 1)
    print('\n\n'+'='*50 + '\n\n LMT C++ class BENCHMARK', Bench)

    Bench = timeit.timeit(setup = setup,stmt = GLMTClass, number = 1)
    print('\n\n'+'='*50 + '\n\n GLMT C++ BENCHMARK', Bench)


def Correctness():
    Phi = np.linspace(-np.pi/2,np.pi/2,100);
    Theta = np.linspace(-np.pi,np.pi,100)

    THETA, PHI = np.meshgrid(Theta, Phi)

    beam = PlaneWave(Wavelength=1e-6)
    BSC = beam.GetBSC(MaxOrder=10)

    argsPy = {"Index":        1.4,
              "Diameter":     1e-6,
              "Wavelength":   1e-6,
              "nMedium":      1.0,
              "Polarization": 0.0,
              "E0":           1.0,
              "R":            1.0,
              "Phi":          Phi,
              "Theta":        Theta}

    argsLMT = {"Index":        1.4,
               "Diameter":     1e-6,
               "Wavelength":   1e-6,
               "nMedium":      1.0,
               "Polarization": 0.0,
               "E0":           1.0}


    argsGLMT = {"Index":        1.4,
                "Diameter":     1e-6,
                "Wavelength":   1e-6,
                "nMedium":      1.0,
                "Polarization": 0.0,
                "E0":           1.0,
                "BSC":          beam._BSC_}

    PyParallel, PyPerpendicular = PyField(**argsPy);

    CppParallel, CppPerpendicular = SPHERELMT(**argsLMT).sFields(Phi=Phi, Theta=Theta, R=1. );

    GLMTCppParallel, GLMTCppPerpendicular = SPHEREGLMT(**argsGLMT).sFields(Phi=Phi, Theta=Theta, R=1. );

    PlotField(Theta, Phi, CppParallel, CppPerpendicular)

    PlotField(Theta, Phi, PyParallel, PyPerpendicular)

    plt.show()



setup = """
import numpy as np
from PyMieSim.LMT.python.Sphere import Fields as PyField

from PyMieSim.Source import PlaneWave
from PyMieSim.LMT.Scatterer import SPHERE as SPHERELMT
from PyMieSim.GLMT.Scatterer import SPHERE as SPHEREGLMT

beam = PlaneWave(Wavelength=1e-6)
BSC = beam.GetBSC(MaxOrder=10)
Phi = np.linspace(-np.pi/2,np.pi/2,800); Theta = np.linspace(-np.pi,np.pi,800)

argsPy = {"Index":        1.4,
          "Diameter":     1e-6,
          "Wavelength":   1e-6,
          "nMedium":      1.0,
          "Polarization": 0.0,
          "E0":           1.0,
          "R":            1.0,
          "Phi":          Phi,
          "Theta":        Theta}

argsLMT = {"Index":        1.4,
           "Diameter":     1e-6,
           "Wavelength":   1e-6,
           "nMedium":      1.0,
           "Polarization": 0.0,
           "E0":           1.0}


argsGLMT = {"Index":        1.4,
            "Diameter":     1e-6,
            "Wavelength":   1e-6,
            "nMedium":      1.0,
            "Polarization": 0.0,
            "E0":           1.0,
            "BSC":          beam._BSC_}

A=SPHERELMT(**argsLMT);
B=SPHEREGLMT(**argsGLMT);

"""

if __name__=='__main__':
    Speed(setup)
    Correctness()











    # -
