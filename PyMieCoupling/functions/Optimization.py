from tqdm import tqdm
import numpy as np
import cupy as cp

from typing import Union
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Operation as Op




def LoopRIDiameter(RIList:       list,
                   DiameterList: list,
                   Detector:     Union[LPmode, Photodiode],
                   QuietMode:    bool = False,
                   cuda:         bool = False,
                   **SKwargs):

    temp = Op.empty(cuda)( [ len(RIList), len(DiameterList) ] )

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes,
                               cuda        = cuda
                               )

            Coupling = Detector.Coupling(Source)

            temp[nr, nd] = Coupling['Parallel']

    return Array( temp )



def ComputeSTD(RIList:       list,
               DiameterList: list,
               Detector:     Union[LPmode, Photodiode],
               QuietMode:    bool = False,
               cuda:         bool = False,
               **SKwargs):

    Array = Array( LoopRIDiameter(RIList,
                                  DiameterList,
                                  Detector,
                                  QuietMode,
                                  cuda,
                                  **SKwargs) )

    return Array



def ComputeMonotonic(RIList:       list,
                     DiameterList: list,
                     Detector:     Union[LPmode, Photodiode],
                     QuietMode:    bool = False,
                     cuda:         bool = False,
                     **SKwargs):

    Array = Array( LoopRIDiameter(RIList,
                                  DiameterList,
                                  Detector,
                                  QuietMode,
                                  cuda,
                                  **SKwargs) )

    val = Monotonic(Array[0,:])


    return val






def Monotonic(Array):
    MaxVal = Array[0]
    N = 0

    for ni, val in enumerate(Array):

        if ni < 1: continue

        if ni == len(Array)-1: continue

        if Array[ni-1] < Array[ni] > Array[ni+1]: N += 1

        if Array[ni-1] > Array[ni] < Array[ni+1]: N += 1


    return N



def Monotonic(Array):

    Grad = np.gradient(Array)

    STD = Grad.std()

    return STD.sum()




class Array(object):
    def __init__(self, Val):
        self.Val = np.array( Val )
        print(self.Val.shape)


    def std(self, arg):
        if arg == 'RI':
            return np.std( self.Val, axis = 0 )

        elif arg == 'Diameter':
            return np.std( self.Val, axis = 1 )

        else:
            Exception('Warning arg put to std() method invalid. Options are ["RI", "Diameter"]')


    def Monotonic(self):

        Grad = np.gradient(self.Val, axis = 1)

        STD = Grad.std( axis = 1)

        return STD[0]



# -
