import importlib
import os
import examples
import logging
import matplotlib._pylab_helpers
import matplotlib.pyplot as plt
from pathlib              import Path
from mayavi               import mlab
from pyface.api           import GUI

from PyMieSim.Directories import StaticPath

def saveMatplotlib():
    """Close the scene."""
    path = os.path.join(StaticPath, f"{name}.png")
    plt.savefig(path)
    plt.close()

def saveMlab():
    """Close the scene."""
    path = os.path.join(StaticPath, f"{name}.png")
    mlab.savefig(path)
    mlab.close()


def SaveCloseMatplotlib():
    """Close matplotlib scene."""
    path = os.path.join(StaticPath, f"{name}.png")

    figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]

    if len(figures) > 0:
        print(f'{name} : Matplotlib example plotting and saving...')
        figures[0].savefig(path)
        plt.close('all')


def SaveCloseMlab():
    """Close mayavi scene."""
    path = os.path.join(StaticPath, f"{name}.png")

    engine = mlab.get_engine()
    scene = engine.current_scene

    if scene is not None:
        print(f'{name} : Mlab example plotting and saving...')
        mlab.savefig(path)
        mlab.close()


def SaveClose():
    """Close all scene."""
    SaveCloseMatplotlib()
    SaveCloseMlab()


def runScript(script):
    global name
    name = script

    try:
        module = importlib.import_module(script)
        module.run(Plot=False, Save=True)

    except:
        logging.warning(f'Script {name} did not conclude, continuing...')

"""
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
"""
runScript('Extra:New-Material-BK7')

runScript('Extra:New-Material-Silver')
