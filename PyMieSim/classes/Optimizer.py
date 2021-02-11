import numpy as np




class Simulator:
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

    def simulate(self, x):
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
            \t Result: {3:.5e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       result)

        if len(self.ParameterName) == 2:
            text = """ \
            Call Number : {0} \
            \t {1}: {2:.5e} \
            \t {3}: {4:.5e} \
            \t Result: {5:.5e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       self.ParameterName[1],
                       x[1],
                       result)

        print(text)
        return result

    def callback(self, xk, *_):
        """Callback function that can be used by optimizers of scipy.optimize.
        The third argument "*_" makes sure that it still works when the
        optimizer calls the callback function with more than one argument. Pass
        to optimizer without arguments or parentheses."""
        s1 = ""
        xk = np.atleast_1d(xk)
        # search backwards in input list for input corresponding to xk
        for i, x in reversed(list(enumerate(self.list_calls_inp))):
            x = np.atleast_1d(x)
            if np.allclose(x, xk):
                break

        for comp in xk:
            s1 += f"{comp:10.5e}\t"
        s1 += f"{self.list_calls_res[i]:10.5e}"

        self.list_callback_inp.append(xk)
        self.list_callback_res.append(self.list_calls_res[i])

        if not self.callback_count:
            s0 = ""
            for j, _ in enumerate(xk):
                tmp = f"Comp-{j+1}"
                s0 += f"{tmp:10s}\t"
            s0 += "Objective"
            print(s0)
        print(s1)
        self.callback_count += 1



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
