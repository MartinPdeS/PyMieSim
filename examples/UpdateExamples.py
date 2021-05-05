import importlib
import examples
from mayavi     import mlab
from pyface.api import GUI
import matplotlib._pylab_helpers
import matplotlib.pyplot as plt

def saveMatplotlib():
    """Close the scene."""
    plt.savefig(f'../docs/images/{name}.png')
    plt.close()

def saveMlab():
    """Close the scene."""
    mlab.savefig(f'../docs/images/{name}.png')
    mlab.close()


def SaveCloseMatplotlib():
    """Close matplotlib scene."""
    figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    if len(figures) >= 1:
        figures[0].savefig(f'../docs/images/{name}.png')
        plt.close('all')


def SaveCloseMlab():
    """Close mayavi scene."""
    engine = mlab.get_engine()
    scene = engine.current_scene
    if scene is not None:
        mlab.savefig(f'../docs/images/{name}.png')
        mlab.close()

def SaveClose():
    """Close all scene."""
    SaveCloseMatplotlib()
    SaveCloseMlab()


def runScript(script):
    global name
    name = script
    module = importlib.import_module(script)

    if module.matplotlib:
        print(f'{name} : Matplotlib example plotting and saving...')
        GUI.invoke_after(2000, SaveClose)
        module.run()

    if module.mlab:
        print(f'{name} : Mlab example plotting and saving...')
        GUI.invoke_after(2000, SaveClose)
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
