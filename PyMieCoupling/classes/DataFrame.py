import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
from tqdm import tqdm

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer
#import cudf
import numpy as np


def GetDataFrame(**kwargs):

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


    def Plot(self, y, Polarization=None):

        if Polarization == 'Parallel':

            self._plot_Parallel(y)


        elif Polarization == 'Perpendicular':

            self._plot_Perpendicular(y)

        elif Polarization == 'Filtered':

            self._plot_Filtered(y)

        else:
            self._plot_Parallel(y)
            self._plot_Perpendicular(y)
            self._plot_Filtered(y)



    def _plot_Parallel(self, y):

        self.ax = self.xs('Parallel').unstack(1).plot(y       = y,
                                                      grid    = True,
                                                      figsize = (8,3),
                                                      title   = '[{1}: {0}] Parallel signal'.format(y, self.DetectorNane),
                                                      ylabel  = y,
                                                      xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)

        plt.subplots_adjust(right=0.8,)



    def _plot_Perpendicular(self, y):
        self.ax = self.xs('Perpendicular').unstack(1).plot(y       = y,
                                                           grid    = True,
                                                           figsize = (8,3),
                                                           title   = '[{1}: {0}] Perpendicular signal'.format(y, self.DetectorNane),
                                                           ylabel  = y,
                                                           xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)


        plt.subplots_adjust(right=0.8,)




    def _plot_Filtered(self, y):
        self.xs('Filtered').unstack(1).plot(y        = y,
                                             grid    = True,
                                             figsize = (8,3),
                                             title   = '[{1}: {0}] Filtered signal'.format(y, self.DetectorNane),
                                             ylabel  = y,
                                             xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)

        plt.subplots_adjust(right=0.8,)



def Frame(RIList:       list,
          DiameterList: list,
          Detector:     Union[LPmode, Photodiode],
          QuietMode:    bool = False,
          **SKwargs):

    Polarization = ['Parallel', 'Perpendicular', 'Filtered']

    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','RI',])

    df = GetDataFrame(index = MI, columns = ['Coupling'])

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(Diameter    = Diameter,
                               Index       = RI,
                               Source      = SKwargs['Source'],
                               Meshes      = Detector.Meshes)

            Coupling = Detector.Coupling(Source)

            df.at[('Parallel', Diameter, RI),'Coupling'] = Coupling['Parallel']

            df.at[('Perpendicular', Diameter, RI),'Coupling'] = Coupling['Perpendicular']

            df.at[('Filtered', Diameter, RI),'Coupling'] = Coupling['Filtered']

    df.Coupling = df.Coupling.astype(float)

    df['Mean'] = df.groupby(['Polarization','Diameter']).Coupling.transform('mean')

    df['STD'] = df.groupby(['Polarization','Diameter']).Coupling.transform('std')

    df.ParaMax = df.xs('Parallel').Coupling.max()

    df.ParaMin = df.xs('Parallel').Coupling.min()

    df.PerpMax = df.xs('Perpendicular').Coupling.max()

    df.PerpMin = df.xs('Perpendicular').Coupling.min()

    df.FilteredMax = df.xs('Filtered').Coupling.max()

    df.FilteredMin = df.xs('Filtered').Coupling.min()

    df.ParaDiff = (df.ParaMax - df.ParaMin).__abs__()

    df.PerpDiff = (df.PerpMax - df.PerpMin).__abs__()

    df.FilteredDiff = (df.FilteredMax - df.FilteredMin).__abs__()

    df.DetectorNane = Detector._name

    return df





# -
