import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"



class ExperimentalDataFrame(pd.DataFrame):

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


    def Plot(self, y='Coupling', **kwargs):

        ax = self.unstack(level=[-3,-1]).plot(y       = y,
                                              grid    = True,
                                              figsize = (8,4),
                                              xlabel  = r'Scatterer diameter [m]',
                                              ylabel  = r'Coupling [u.a.]',
                                              **kwargs)

        ax.legend(prop={'size': 8})
        plt.show()


class S1S2DataFrame(pd.DataFrame):

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


    def Plot(self, **kwargs):

        ax = self.unstack(level=[0,1]).plot(y       = 'S1',
                                            grid    = True,
                                            figsize = (8,4),
                                            xlabel  = r'$\phi$ angle [degree]',
                                            ylabel  = r'$|S1|$',
                                            **kwargs)

        ax1 = self.unstack(level=[0,1]).plot(y       = 'S2',
                                             grid    = True,
                                             figsize = (8,4),
                                             xlabel  = r'$\phi$ angle [degree]',
                                             ylabel  = r'$|S2|$',
                                             **kwargs)

        ax.legend(prop={'size': 8})

        ax1.legend(prop={'size': 8})

        plt.show()




class QscaDataFrame(pd.DataFrame):

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


    def Plot(self, **kwargs):

        ax = self.unstack(level=[1]).plot(y       = 'Qsca',
                                          grid    = True,
                                          figsize = (8,4),
                                          xlabel  = r'Scatterer diameter [m]',
                                          ylabel  = r'Q$_{Scattering}$',
                                          **kwargs)

        ax.legend(prop={'size': 8})


        plt.show()
















# -
