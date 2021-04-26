import os

import PyMieSim

RootPath = PyMieSim.__path__[0]

MaterialPath = os.path.join(RootPath, 'Data/_Material')

NPZPath = os.path.join(MaterialPath, 'npz')

LPModePath = os.path.join(RootPath, 'Data/LPmodes')
