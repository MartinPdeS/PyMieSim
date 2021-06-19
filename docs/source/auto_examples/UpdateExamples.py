"""
Update Examples
"""

import importlib
import os
import examples
import logging
import matplotlib._pylab_helpers
import matplotlib.pyplot as plt
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

    if len(figures) >= 1:
        figures[0].savefig(path)
        plt.close('all')


def SaveCloseMlab():
    """Close mayavi scene."""
    path = os.path.join(StaticPath, f"{name}.png")

    engine = mlab.get_engine()
    scene = engine.current_scene
    if scene is not None:
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

        #if module.matplotlib:
        print(f'{name} : Matplotlib example plotting and saving...')
        GUI.invoke_after(200, SaveClose)
        module.run()

        #if module.mlab:
        print(f'{name} : Mlab example plotting and saving...')
        GUI.invoke_after(200, SaveClose)
        module.run()
    except:
        logging.warning(f'Script {name} did not conclude, continuing...')


runScript('LPMode')

runScript('Photodiode')

runScript('S1S2')

runScript('Stokes')

runScript('SPF')

runScript('FarField')

runScript('Mie-resonances')

runScript('Qscattering')

runScript('Coupling-vs-wavelength')

runScript('Qsca-vs-diameter')

runScript('Coupling-vs-diameter')

runScript('Experiment:Goniometer')

runScript('New-Material')

runScript('New-Material-Silver')
