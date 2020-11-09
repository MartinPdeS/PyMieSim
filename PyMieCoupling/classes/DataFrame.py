import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
from tqdm import tqdm

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Scattering import Scatterer
#import cudf


def GetDataFrame(cuda, **kwargs):
    if cuda:
        return DataFrameCPU(**kwargs)
    else:
        return DataFrameCPU(**kwargs)


class DataFrameCPU(pd.DataFrame):

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
                                            xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)

        plt.subplots_adjust(right=0.8,)

        self.xs('Perpendicular').unstack(1).plot(y       = kwargs['y'],
                                                 grid    = True,
                                                 figsize = (8,3),
                                                 title   = '[{1}: {0}] Perpendicular field'.format(kwargs['y'], self.DetectorNane),
                                                 ylabel  = 'Coupling',
                                                 xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)

        plt.subplots_adjust(right=0.8,)





class DataFrameGPU(pd.DataFrame):

    def __init__(self,**kwargs):
        cudf.DataFrame.__init__(self,**kwargs)

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
                                            xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)
        plt.subplots_adjust(right=0.8,)

        self.xs('Perpendicular').unstack(1).plot(y       = kwargs['y'],
                                                 grid    = True,
                                                 figsize = (8,3),
                                                 title   = '[{1}: {0}] Perpendicular field'.format(kwargs['y'], self.DetectorNane),
                                                 ylabel  = 'Coupling',
                                                 xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)
        plt.subplots_adjust(right=0.8,)





def Frame(RIList:       list,
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





# -
