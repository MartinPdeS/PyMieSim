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

def runScript(script):
    global name
    name = script
    module = importlib.import_module(script)

    if module.matplotlib:
        print('Matplotlib example plotting and saving...')
        GUI.invoke_after(2000, saveMatplotlib)
        module.run()

    if module.mlab:
        print('Mlab example plotting and saving...')
        GUI.invoke_after(2000, saveMlab)
        module.run()


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

runScript('New-Material')

runScript('New-Material-Silver')
