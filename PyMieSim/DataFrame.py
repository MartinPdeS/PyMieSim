import pandas as pd
import matplotlib.pyplot as plt

def show(func):
    def inner(*args, **kwargs):
        try: return func(*args, **kwargs)

        finally: plt.show()

    return inner



class ExperimentalDataFrame(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self, **kwargs)
        self.ax = None

    @property
    def Parallel(self):
        return self.xs('Parallel')

    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')

    @show
    def Plot(self, y='Coupling', **kwargs):

        fig = self.unstack(level=[-3,-1]).plot(y       = y,
                                              grid    = True,
                                              figsize = (8,4),
                                              xlabel  = r'Scatterer diameter [m]',
                                              ylabel  = r'Coupling [Watt]',
                                              **kwargs)

        fig.legend(prop={'size': 8}, title='Detec. , Refrac. index',loc=2)

        return fig



class S1S2DataFrame(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)
        self.ax = None

    @property
    def Parallel(self):
        return self.xs('Parallel')

    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')

    @show
    def Plot(self, **kwargs):

        fig = self.unstack(level=[0,1]).plot(y       = 'S1',
                                            grid    = True,
                                            figsize = (8,4),
                                            xlabel  = r'$\phi$ angle [degree]',
                                            ylabel  = r'$|S1|$',
                                            **kwargs)

        fig1 = self.unstack(level=[0,1]).plot(y       = 'S2',
                                             grid    = True,
                                             figsize = (8,4),
                                             xlabel  = r'$\phi$ angle [degree]',
                                             ylabel  = r'$|S2|$',
                                             **kwargs)

        fig.legend(prop={'size': 8})

        fig1.legend(prop={'size': 8})

        return(fig, fig1)



class QscaDataFrame(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)
        self.ax = None

    @property
    def Parallel(self):
        return self.xs('Parallel')

    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')

    @show
    def Plot(self, **kwargs):

        fig = self.unstack(level=[1]).plot(y       = 'Qsca',
                                          grid    = True,
                                          figsize = (8,4),
                                          xlabel  = r'Scatterer diameter [m]',
                                          ylabel  = r'Q$_{Scat}$ [Watt.m$^{-2}$]',
                                          **kwargs)

        fig.legend(prop={'size': 8}, title='Refrac. index', loc=4)

        return fig

















# -
