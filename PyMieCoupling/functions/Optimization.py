from tqdm import tqdm
import numpy as np

from typing import Union
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Scattering import Scatterer




def LoopRIDiameter(RIList:       list,
                   DiameterList: list,
                   Detector:     Union[LPmode, Photodiode],
                   QuietMode:    bool = False,
                   Polarization: str  = 'Parallel',
                   **SKwargs):

    temp = np.empty( [ len(RIList), len(DiameterList) ] )


    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes
                               )

            Coupling = Detector.Coupling(Source)

            temp[nr, nd] = Coupling[Polarization]

    return Array(temp)



class Array(np.ndarray):
    def __new__(cls, *args, **kwargs):
        this = np.array(*args, **kwargs)
        this = np.asarray(this).view(cls)
        return this

    def __array_finalize__(self, obj):
        pass


    def __init__(self, arr):
        self = np.array([1,2,3])


    def Cost(self, arg='RI'):
        if arg == 'RI':
            return self.std(axis=0).sum()

        if arg == 'RI/Mean':
            return self.std(axis=0).sum()/self.mean()

        if arg == 'Monotonic':
            return self.Monotonic()

        if arg == 'Mean':
            return -self.mean()

        if arg == 'Max':
            return -self.max()

        if arg == 'Min':
            return self.max()


    def Monotonic(self):

        Grad = np.gradient(self, axis = 1)

        STD = Grad.std( axis = 1)

        return STD[0]



# -
