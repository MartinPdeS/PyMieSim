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
            if   'max' in arg:  return np.max(self)
            elif 'min' in arg:  return np.min(self)
            elif 'mean' in arg: return np.mean(self)

        if len(arg) == 2:
            if   arg[0] == 'rsd':        func = self.rsd
            elif arg[0] == 'monotonic':  func = self.Monotonic

            if   arg[1] == 'ri':           return np.mean( func(self, axis = 4) )
            elif arg[1] == 'diameter':     return np.mean( func(self, axis = 3) )
            elif arg[1] == 'polarization': return np.mean( func(self, axis = 2) )
            elif arg[1] == 'wavelength':   return np.mean( func(self, axis = 1) )
            elif arg[1] == 'detector':     return np.mean( func(self, axis = 0) )

        raise ValueError(f"Invalid metric input. \nList of metrics: {MetricList}")



    def Monotonic(self, axis):

        axis = axis.lower()

        arr = np.gradient(self.data,
                          axis = self.conf['order'][axis]).std( axis = self.conf['order'][axis])

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Mean(self, axis):

        axis = axis.lower()

        arr = np.mean(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Std(self, axis):

        axis = axis.lower()

        arr = np.std(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def Rsd(self, axis):

        axis = axis.lower()

        arr = np.std(self.data, axis=self.conf['order'][axis] ) \
             /np.mean(self.data, axis=self.conf['order'][axis] )

        conf = self.UpdateConf(axis)

        return PMSArray(array=arr, conf=conf)


    def UpdateConf(self, axis):

        newConf = copy.deepcopy(self.conf)

        newConf['order'].pop(axis)
        newConf['dimension'].pop(axis)

        for n, key in enumerate(newConf['order'].keys()):
            newConf['order'][key] = n

        return newConf


    def Plot(self, x):
        fig = plt.figure(figsize=(8,4))
        x = x.lower()
        shape = list(self.data.shape)
        for key, order in self.conf['order'].items():
            if x == key:
                shape[order] = None
                xlabel = self.conf['label'][key]
                xval   = self.conf['dimension'][key]


        for idx in product(*[range(s) if s is not None else [slice(None)] for s in shape]):
            plt.plot(xval,
                     self.data[idx],
                     label = self.GetLabel(x, idx))

        plt.xlabel(xlabel)
        plt.ylabel(self.conf['label']['variable'])
        plt.grid()
        plt.legend(fontsize=8)
        plt.show()


    def GetLabel(self, x, idx):
        label = ''

        for key in self.conf['order']:

            if x != key:
                index = idx[self.conf['order'][key]]
                val = self.conf['dimension'][key][index]
                label += f"{key[:3]}.:{val} | "

        return label


    def __getitem__(self, key):
        return self.data[key]


    def __setitem__(self, key, value):
        self.data[key] = value


    def __str__(self):
        print(self.conf['name'])
        text = f'PyMieArray \nVariable: {self.Name}\n' + '='*90 + '\n'
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

        self.dim = { 'detector'      : True,
                      'wavelength'   : True,
                      'polarization' : True,
                      'diameter'     : True,
                      'index'        : True}


    def Cost(self, arg = 'max'):

        arg = arg.lower().split('+', 2)

        if len(arg) == 1:
            if   'max' in arg:  return np.max(self)
            elif 'min' in arg:  return np.min(self)
            elif 'mean' in arg: return np.mean(self)

        if len(arg) == 2:
            if   arg[0] == 'rsd':        func = self.rsd
            elif arg[0] == 'monotonic':  func = self.Monotonic

            if   arg[1] == 'ri':           return np.mean( func(self, axis = 4) )
            elif arg[1] == 'diameter':     return np.mean( func(self, axis = 3) )
            elif arg[1] == 'polarization': return np.mean( func(self, axis = 2) )
            elif arg[1] == 'wavelength':   return np.mean( func(self, axis = 1) )
            elif arg[1] == 'detector':     return np.mean( func(self, axis = 0) )

        raise ValueError(f"Invalid metric input. \nList of metrics: {MetricList}")


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
