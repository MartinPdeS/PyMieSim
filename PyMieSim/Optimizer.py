import numpy as np
from scipy.optimize import minimize
from mayavi import mlab



class Optimize:
    def __init__(self,
                 ExperimentalSet,
                 Metric,
                 Parameter,
                 X0,
                 WhichDetector,
                 MinVal,
                 MaxVal,
                 FirstStride,
                 MaxIter=50,
                 Tol=1e-10):

        self.ExperimentalSet = ExperimentalSet
        self.Metric          = Metric
        self.Parameters       = Parameter
        self.X0              = X0
        self.WhichDetector   = WhichDetector
        self.MinVal          = MinVal
        self.MaxVal          = MaxVal
        self.FirstStride     = FirstStride
        self.MaxIter         = MaxIter
        self.Tol             = Tol

        if len(self.Parameters)   == 1: self.Result = self.Optimize1()

        elif len(self.Parameters) == 2: self.Result = self.Optimize2()


    def Optimize1(self):

        def EvalFunc(x):
            Penalty = 0
            for n in range(len(self.Parameters)):
                if self.MinVal[n] and x[0]< self.MinVal[n]: Penalty += np.abs( x[0]*100 ); x[0] = self.MinVal[n]
                if self.MinVal[n] and x[0]> self.MaxVal[n]: Penalty += np.abs( x[0]*100 ); x[0] = self.MaxVal[n]

            setattr(self.ExperimentalSet.Detectors[self.WhichDetector], self.Parameters[0], x[0])

            return self.ExperimentalSet.Coupling.Cost(self.Metric) + Penalty

        Minimizer = Caller(EvalFunc, ParameterName = self.Parameters)

        return minimize(fun      = Minimizer.optimize,
                        x0       = self.X0,
                        method   = 'COBYLA',
                        tol      = self.Tol,
                        options  = {'maxiter': self.MaxIter, 'rhobeg':self.FirstStride})


    def Optimize2(self):

        def EvalFunc(x):
            Penalty = 0
            for n in range(len(self.Parameters)):
                if self.MinVal[n] and x[0]< self.MinVal[n]: Penalty += np.abs( x[0]*100 ); x[0] = self.MinVal[n]
                if self.MinVal[n] and x[0]> self.MaxVal[n]: Penalty += np.abs( x[0]*100 ); x[0] = self.MaxVal[n]

            setattr(self.ExperimentalSet.Detectors[self.WhichDetector], self.Parameters[0], x[0])
            setattr(self.ExperimentalSet.Detectors[self.WhichDetector], self.Parameters[1], x[1])

            return self.ExperimentalSet._Coupling(self.WhichDetector).Cost(self.Metric) + Penalty

        Minimizer = Caller(EvalFunc, ParameterName = self.Parameters)

        return minimize(fun      = Minimizer.optimize,
                        x0       = self.X0,
                        method   = 'COBYLA',
                        tol      = self.Tol,
                        options  = {'maxiter': self.MaxIter, 'rhobeg':self.FirstStride})


class Caller:
    def __init__(self, function, ParameterName: list):
        self.ParameterName = ParameterName
        self.f = function # actual objective function
        self.num_calls = 0 # how many times f has been called
        self.callback_count = 0 # number of times callback has been called, also measures iteration count
        self.list_calls_inp = [] # input of all calls
        self.list_calls_res = [] # result of all calls
        self.decreasing_list_calls_inp = [] # input of calls that resulted in decrease
        self.decreasing_list_calls_res = [] # result of calls that resulted in decrease
        self.list_callback_inp = [] # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = [] # only appends results on callback, as such they correspond to the iterations

    def optimize(self, x):
        """Executes the actual simulation and returns the result, while
        updating the lists too. Pass to optimizer without arguments or
        parentheses."""
        result = self.f(x) # the actual evaluation of the function
        if not self.num_calls: # first call is stored in all lists
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
            self.list_callback_inp.append(x)
            self.list_callback_res.append(result)
        elif result < self.decreasing_list_calls_res[-1]:
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
        self.list_calls_inp.append(x)
        self.list_calls_res.append(result)
        self.num_calls += 1


        if len(self.ParameterName) == 1:

            text = """ \
            Call Number : {0} \
            \t {1}: {2:.5e} \
            \t Result: {3:.10e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       result)

        if len(self.ParameterName) == 2:
            text = """ \
            Call Number : {0} \
            \t {1}: {2:.5e} \
            \t {3}: {4:.5e} \
            \t Result: {5:.10e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       self.ParameterName[1],
                       x[1],
                       result)

        print(text)
        return result




class OptArray(np.ndarray):
    def __new__(cls, *args, **kwargs):
        this = np.array(*args, **kwargs, copy=False)
        this = np.asarray(this).view(cls)
        return this

    def __array_finalize__(self, obj):
        pass


    def __init__(self, arr):
        pass


    def Cost(self, arg='RI'):
        if arg == 'RI_STD':
            return self.std(axis=0).sum()

        if arg == 'RI_RSD':
            return self.std(axis=0).sum()/self.mean()

        if arg == 'Size_STD':
            return self.std(axis=1).sum()

        if arg == 'Size_RSD':
            return self.std(axis=1).sum()/self.mean()

        if arg == 'Monotonic':
            return self.Monotonic()

        if arg == 'Mean':
            return -self.mean()

        if arg == 'Max':
            return -self.max()

        if arg == 'Min':
            return self.max()


    def Monotonic(self):

        Grad = np.gradient(self, axis = 1)

        STD = Grad.std( axis = 1)

        return STD[0]
