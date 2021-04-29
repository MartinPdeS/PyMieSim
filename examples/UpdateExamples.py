import importlib
import examples
from mayavi     import mlab
from pyface.api import GUI
import matplotlib.pyplot as plt

def saveMatplotlib():
    """Close the scene."""
    plt.savefig(f'docs/images/{name}.png')
    plt.close()

def saveMlab():
    """Close the scene."""
    mlab.savefig(f'docs/images/{name}.png')
    mlab.close()

def runScriptMlab(script):
    global name
    name = script
    GUI.invoke_after(2000, saveMlab)
    module = importlib.import_module(script)

def runScriptMatplotlib(script):
    global name
    name = script
    GUI.invoke_after(2000, saveMatplotlib)
    module = importlib.import_module(script)


runScriptMlab('LPMode')

runScriptMlab('Photodiode')

runScriptMatplotlib('S1S2')

runScriptMlab('Stokes')

runScriptMlab('SPF')

runScriptMlab('FarField')

runScriptMatplotlib('Mie-resonances')

runScriptMatplotlib('Qscattering')

runScriptMatplotlib('Coupling-vs-wavelength')

runScriptMatplotlib('Coupling-vs-diameter')

runScriptMatplotlib('New-Material')

runScriptMatplotlib('New-Material-Silver')
