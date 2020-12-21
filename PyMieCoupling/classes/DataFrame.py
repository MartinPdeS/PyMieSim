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


    def Plot(self, y, Scale='Linear'):

        for Polar in self.attrs['Filter']:
            self._plot(y, Polar, Scale)


    def _plot(self, y, Filter, Scale):

        self.ax = self.xs(Filter).unstack(1).plot(y       = y,
                                                  grid    = True,
                                                  figsize = (8,3.5),
                                                  title   = r'[{0}] Filter: {1} [Degree]'.format(self.DetectorNane, Filter),
                                                  ylabel  = y,
                                                  xlabel  = r'Scatterer diameter [m]')

        self.ax.tick_params(labelsize='small')
        self.ax.legend(bbox_to_anchor=(1, 1), ncol=1)

        if Scale == 'Logarithmic':
            self.ax.set_yscale('log')

        plt.subplots_adjust(right=0.8,)

        plt.show(block=False)
