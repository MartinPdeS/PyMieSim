import pandas as pd
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





# -
