import numpy as np
import copy
import matplotlib.pyplot as plt
from itertools import product

from PyMieSim.Config import MetricList



class PMSArray(object):

    def __init__(self, array, conf):
        self.data = array
        self.conf = conf


    def Cost(self, arg = 'max'):

        arg = arg.lower().split('+', 2)

        if len(arg) == 1:
            if   'max'  in arg : return np.max(self)
            elif 'min'  in arg : return np.min(self)
            elif 'mean' in arg : return np.mean(self)

        if len(arg) == 2:
            if   arg[0] == 'rsd'          : func = self.rsd
            elif arg[0] == 'monotonic'    : func = self.Monotonic

            if   arg[1] == 'ri'           : return np.mean( func(self.data, axis = 4) )
            elif arg[1] == 'diameter'     : return np.mean( func(self.data, axis = 3) )
            elif arg[1] == 'polarization' : return np.mean( func(self.data, axis = 2) )
            elif arg[1] == 'wavelength'   : return np.mean( func(self.data, axis = 1) )
            elif arg[1] == 'detector'     : return np.mean( func(self.data, axis = 0) )

        raise ValueError(f"Invalid metric input. \nList of metrics: {MetricList}")



    def Monotonic(self, axis):
        """Method compute and the monotonic value of specified axis.
        The method then return a new PMSArray daughter object compressed in
        the said axis.

        Parameters
        ----------
        axis : :class:`str`
            Axis for which to perform the operation.

        Returns
        -------
        :class:`PMSArray`
            New PMSArray instance containing the monotonic metric value of axis.

        """
        axis = axis.lower()

        arr  = np.gradient(self.data,
                           axis = self.conf['order'][axis])\
                           .std( axis = self.conf['order'][axis])

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Mean(self, axis):
        """Method compute and the mean value of specified axis.
        The method then return a new PMSArray daughter object compressed in
        the said axis.

        Parameters
        ----------
        axis : :class:`str`
            Axis for which to perform the operation.

        Returns
        -------
        :class:`PMSArray`
            New PMSArray instance containing the mean value of axis.

        """
        axis = axis.lower()

        arr  = np.mean(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Std(self, axis):
        """Method compute and the std value of specified axis.
        The method then return a new PMSArray daughter object compressed in
        the said axis.

        Parameters
        ----------
        axis : :class:`str`
            Axis for which to perform the operation.

        Returns
        -------
        :class:`PMSArray`
            New PMSArray instance containing the std value of axis.

        """
        axis = axis.lower()

        arr  = np.std(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Rsd(self, axis):
        """Method compute and the rsd value of specified axis.
        The method then return a new PMSArray daughter object compressed in
        the said axis.
        rsd is defined as std/mean.

        Parameters
        ----------
        axis : :class:`str`
            Axis for which to perform the operation.

        Returns
        -------
        :class:`PMSArray`
            New PMSArray instance containing the rsd value of axis.

        """
        axis = axis.lower()

        arr  = np.std(self.data, axis=self.conf['order'][axis] ) \
              /np.mean(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def UpdateConf(self, axis):

        newConf = self.conf.copy()#copy.copy(self.conf)

        newConf['order'].pop(axis)
        newConf['dimension'].pop(axis)

        for n, key in enumerate(newConf['order'].keys()):
            newConf['order'][key] = n

        return newConf


    def Plot(self, Scale = 'linear', *args, **kwargs):
        if Scale.lower() in ['lin', 'linear']      : plot = plt.plot
        if Scale.lower() in ['log', 'logarithmic'] : plot = plt.loglog;

        if self.conf['variable']['name'] == 'Efficiencies':
            return self.PlotEfficiencies(*args, **kwargs, plot=plot)

        if self.conf['variable']['name'] == 'Coupling':
            return self.PlotCoupling(*args, **kwargs, plot=plot)



    def PlotCoupling(self, x, plot, Testing=False, *args, **kwargs):

        fig   = self.PrepareFigure()
        x     = x.lower()

        DimSlicer, xlabel, xval = self.GetSlicer(x)

        for idx in product(*DimSlicer):
            plot(xval,
                 self.data[idx],
                 label = self.GetLegend(x, idx),
                 *args,
                 **kwargs)

        plt.xlabel(xlabel)
        plt.legend(fontsize=8)

        if Testing == False:
            plt.show()
        else:
            plt.close('all')


    def GetSlicer(self, x):
        shape = list(self.data.shape)

        for key, order in self.conf['order'].items():
            if key.lower() == x.lower():
                shape[order] = None
                xlabel       = self.conf['label'][key]
                xval         = self.conf['dimension'][key]

        DimSlicer = [range(s) if s is not None else [slice(None)] for s in shape]

        return DimSlicer, xlabel, xval


    def PrepareFigure(self):
        fig   = plt.figure(figsize=(8,4))

        plt.grid()

        plt.ylabel(self.conf['variable']['name'] + ' ' + self.conf['variable']['unit'])

        return fig


    def PlotEfficiencies(self, x, plot, Testing=False,  *args, **kwargs):

        fig   = self.PrepareFigure()

        x     = x.lower()

        DimSlicer, xlabel, xval = self.GetSlicer(x)

        for ni, idx in enumerate( product( *DimSlicer ) ):
            plot(xval,
                 self.data[idx],
                 label = self.conf['variable']['namelist'][idx[-1]] + ' | ' + self.GetLegend(x, idx),
                 *args,
                 **kwargs)

        plt.xlabel(xlabel)

        plt.legend(fontsize=8)

        if Testing == False:
            plt.show()
        else:
            plt.close('all')


    def GetLegend(self, axis, idx):
        """Method generate and return the legend text for the specific plot.

        Parameters
        ----------
        axis : :class:`str`
            Axis which is used for x-axis of the specific plot
        idx : :class:`tuple`
            Dimension indices of the specific plot

        Returns
        -------
        :class:`str`
            Text for the legend

        """
        label = ''

        for key in self.conf['order']:

            if axis.lower() != key.lower():

                if key == 'material':
                    index  = idx[self.conf['order'][key]]
                    val    = self.conf['dimension'][key][index]
                    format = self.conf['format']['material']
                    label += f"{key}: { val } | "

                else:
                    index  = idx[self.conf['order'][key]]
                    val    = self.conf['dimension'][key][index]
                    format = self.conf['format'][key]
                    label += f"{key}= {val:{format}} | "

        return label


    def __getitem__(self, key):
        return self.data[key]


    def __setitem__(self, key, value):
        self.data[key] = value


    def __str__(self):
        name = self.conf['name']
        text =  f'PyMieArray \nVariable: {name}\n' + '='*90 + '\n'
        text += f"{'Parameter':13s}\n" + '-'*90 + '\n'
        for key, val in self.conf['order'].items():
            text += f"""{key:13s}\
                        | dimension = {val:2d}\
                        | size = {len(self.conf['dimension'][key]):2d}\
                         \n"""

        text += '='*90 + '\n'
        return text


class Opt5DArray(np.ndarray):
    def __new__(cls, *args, **kwargs):
        this = np.array(*args, **kwargs, copy=False)
        this = np.asarray(this).view(cls)

        return this


    def __array_finalize__(self, obj):
        pass


    def __init__(self, arr, Name=''):
        self.Name         = Name

        self.dim = { 'detector'     : True,
                     'wavelength'   : True,
                     'polarization' : True,
                     'diameter'     : True,
                     'index'        : True}

    def DefineCostFunc(self, arg):
        arg = arg.lower().split('+', 2)

        if len(arg) == 1:
            if   'max'  in arg : self.CostFunc = np.max
            elif 'min'  in arg : self.CostFunc = np.min
            elif 'mean' in arg : self.CostFunc = np.mean

        if len(arg) == 2:
            if   arg[0] == 'rsd'       : func = self.rsd
            elif arg[0] == 'monotonic' : func = self.Monotonic
            elif arg[0] == 'max'       : func = np.max
            elif arg[0] == 'min'       : func = np.min
            elif arg[0] == 'mean'      : func = np.mean


            if   arg[1] == 'all'          : self.CostFunc = np.mean( func(self) )
            elif arg[1] == 'ri'           : self.CostFunc = np.mean( func(self, axis = 4) )
            elif arg[1] == 'diameter'     : self.CostFunc = np.mean( func(self, axis = 3) )
            elif arg[1] == 'polarization' : self.CostFunc = np.mean( func(self, axis = 2) )
            elif arg[1] == 'wavelength'   : self.CostFunc = np.mean( func(self, axis = 1) )
            elif arg[1] == 'detector'     : self.CostFunc = np.mean( func(self, axis = 0) )

            raise ValueError(f"Invalid metric input. \nList of metrics: {MetricList}")



    def Cost(self):

        return self.CostFunc(self)




    def Monotonic(self, axis):

        Grad = np.gradient(self, axis = axis)

        STD = Grad.std( axis = axis)

        return STD[0]


    def rsd(self, array, axis):
        return np.std(array, axis)/np.mean(array, axis)


    def RIMonotonic(self):

        Grad = np.gradient(self, axis = 0)

        STD = Grad.std( axis = 0)

        return STD[0]
