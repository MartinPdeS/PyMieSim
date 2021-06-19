
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import importlib
import os
import examples
import logging
import matplotlib._pylab_helpers
import time
import matplotlib.pyplot as plt
from pathlib              import Path
from mayavi               import mlab
from pyface.api           import GUI



from PyMieSim.Tools.Directories import StaticPath


def runScript(script):

    dir = f'docs/images/{Path(__file__).stem}'

    dir = os.path.join(StaticPath, dir)

    try:
        module = importlib.import_module(script)
        module.run(Plot=False, Save=True, Directory=dir)

    except:
        logging.warning(f'Script {script} did not conclude, continuing...')


runScript('Coupling:LPMode')

runScript('Coupling:LPMode')

runScript('Coupling:Photodiode')

runScript('Fields:S1S2')

runScript('Fields:Stokes')

runScript('Fields:SPF')

runScript('Fields:FarField')

runScript('Experiment:Mie-resonances')

runScript('Experiment:Qscattering')

runScript('Experiment:Coupling-vs-wavelength')

runScript('Experiment:Qsca-vs-diameter')

runScript('Experiment:Coupling-vs-diameter')

runScript('Experiment:Goniometer')

runScript('Extra:New-Material-BK7')

runScript('Extra:New-Material-Silver')
