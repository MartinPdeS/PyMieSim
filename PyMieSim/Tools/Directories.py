import os

import PyMieSim

RootPath     = PyMieSim.__path__[0]

ZeroPath     = os.path.dirname(RootPath)

MaterialPath = os.path.join(RootPath, 'Data/_Material')

NPZPath      = os.path.join(MaterialPath, 'npz')

StaticPath   = os.path.join(ZeroPath, 'docs/images')

LPModePath   = os.path.join(RootPath, 'Data/LPmodes')

RTDExample   = 'https://pymiesim.readthedocs.io/en/latest/Examples.html'

RTDMaterial  = 'https://pymiesim.readthedocs.io/en/latest/Material.html'

RTDLPMode    = 'https://pymiesim.readthedocs.io/en/latest/LPModes.html'
