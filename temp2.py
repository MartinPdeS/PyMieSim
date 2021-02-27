from PyMieSim.LMT.PyBind.Sphere import PyFieldsNoPolarization, PyFields

a = PyFields(1.4,
             1e-6,
             1e-6,
             1,
             [1,2,3],
             [1,2,3],
             3,0)
print(a)
