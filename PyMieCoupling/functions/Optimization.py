import numpy as np
from tqdm import tqdm
import pandas as pd
from typing import Union
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer


class WrapDataFrame(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)

    @property
    def Parallel(self):
        return self.xs('Parallel')

    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')

    def plot(self, **kwargs):
        self.xs('Parallel').unstack(1).plot(y       = kwargs['y'],
                                            grid    = True,
                                            figsize = (8,3),
                                            title   = '[{1}: {0}] Parallel field'.format(kwargs['y'], self.DetectorNane),
                                            ylabel  = 'Coupling',
                                            #color   = 'C1',
                                            xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)
        plt.subplots_adjust(right=0.8,)

        self.xs('Perpendicular').unstack(1).plot(y       = kwargs['y'],
                                                 grid    = True,
                                                 figsize = (8,3),
                                                 title   = '[{1}: {0}] Perpendicular field'.format(kwargs['y'], self.DetectorNane),
                                                 ylabel  = 'Coupling',
                                                 #color   = 'CO',
                                                 xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)
        plt.subplots_adjust(right=0.8,)


def CouplingStat(RIList: list,
                 DiameterList: list,
                 Detector: Union[LPmode, Photodiode],
                 QuietMode: bool = False,
                 **SKwargs) -> WrapDataFrame:

    Polarization = ['Parallel', 'Perpendicular']

    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','RI',])

    df = WrapDataFrame(index = MI, columns = ['Coupling'])

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes,
                               GPU         = SKwargs['GPU']
                               )

            Perp, Para = PointFieldCoupling(Detector = Detector, Source = Source)

            df.at[('Parallel', Diameter, RI),'Coupling'] = Perp

            df.at[('Perpendicular', Diameter, RI),'Coupling'] = Para

    df.Coupling = df.Coupling.astype(float)

    df['Mean'] = df.groupby(['Polarization','Diameter']).Coupling.transform('mean')

    df['STD'] = df.groupby(['Polarization','Diameter']).Coupling.transform('std')

    df.ParaMax = df.xs('Parallel').Coupling.max()

    df.ParaMin = df.xs('Parallel').Coupling.min()

    df.PerpMax = df.xs('Perpendicular').Coupling.max()

    df.PerpMin = df.xs('Perpendicular').Coupling.min()

    df.ParaDiff = np.abs(df.ParaMax - df.ParaMin)

    df.PerpDiff = np.abs(df.PerpMax - df.PerpMin)

    df.DetectorNane = Detector._name

    return df


def OptimizeMono(DiameterList: list,
                 Detector):
    pass
