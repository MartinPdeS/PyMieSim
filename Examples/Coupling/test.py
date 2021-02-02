from PyMieCoupling.cpp.interface import GetS1S2Qsca as S1S2_CPP
import numpy as np

AngleList = np.linspace(-np.pi/2,np.pi/2,101)#.tolist()

SizeParam, index = 0.2, 1.4


S1S2, Qsca = S1S2_CPP(index, SizeParam, AngleList);

print(Qsca)
