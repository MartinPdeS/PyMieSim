from tqdm import tqdm
import numpy as np
import cupy as cp
import pandas as pd
from typing import Union
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.DataFrame import GetDataFrame


def CouplingStat(RIList:       list,
                 DiameterList: list,
                 Detector:     Union[LPmode, Photodiode],
                 QuietMode:    bool = False,
                 cuda:         bool = False,
                 **SKwargs):

    Polarization = ['Parallel', 'Perpendicular']

    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','RI',])

    df = GetDataFrame(index = MI, columns = ['Coupling'], cuda = cuda)

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes,
                               cuda        = cuda)

            Coupling = Detector.Coupling(Source)

            df.at[('Parallel', Diameter, RI),'Coupling'] = Coupling['Perpendicular']

            df.at[('Perpendicular', Diameter, RI),'Coupling'] = Coupling['Perpendicular']

    df.Coupling = df.Coupling.astype(float)

    df['Mean'] = df.groupby(['Polarization','Diameter']).Coupling.transform('mean')

    df['STD'] = df.groupby(['Polarization','Diameter']).Coupling.transform('std')

    df.ParaMax = df.xs('Parallel').Coupling.max()

    df.ParaMin = df.xs('Parallel').Coupling.min()

    df.PerpMax = df.xs('Perpendicular').Coupling.max()

    df.PerpMin = df.xs('Perpendicular').Coupling.min()

    df.ParaDiff = (df.ParaMax - df.ParaMin).__abs__()

    df.PerpDiff = (df.PerpMax - df.PerpMin).__abs__()

    df.DetectorNane = Detector._name

    return df




def CompileStat(RIList:       list,
                 DiameterList: list,
                 Detector:     Union[LPmode, Photodiode],
                 QuietMode:    bool = False,
                 cuda:         bool = False,
                 **SKwargs):


    Array = np.empty( [ len(RIList), len(DiameterList) ] )
    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes,
                               cuda        = cuda
                               )

            Coupling = Detector.Coupling(Source)

            Array[nr, nd] = Coupling['Parallel']

    return Array















# -
