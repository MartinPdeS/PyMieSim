import pandas as pd
import matplotlib.pyplot as plt


class DataFrameCPU(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)
        self.Filter = None
        self.ax = None


    @property
    def Parallel(self):
        return self.xs('Parallel')


    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')


    def _Plot(self, y, Scale='Linear'):

        for Detector in self.attrs['Detectors'] .keys():
            print(Detector)
        self._plot(y, Detector, Scale)


    def Plot(self, y, Scale='Linear'):

        self.unstack(level=[-3,-1]).plot(y       = 'Coupling',
                                         grid    = True,
                                         figsize = (8,3.5),
                                         xlabel  = r'Scatterer diameter [m]',
                                         ylabel  = r'Coupling [u.a.]')
















# -
