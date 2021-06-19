import json
import os
import csv
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


PATH = os.path.join( Path(__file__).parent, 'Data/_Material' )

from PyMieSim.Tools.utils       import IO
from PyMieSim.Tools.Directories import ZeroPath
from PyMieSim.BaseClasses       import BaseMaterial


class Material(BaseMaterial):
    def __init__(self, name):
        self._Data = None
        self.__name__ = name

        with open(os.path.join(PATH, 'Meta.json'), 'r+' ) as f:
            META                     = json.load(f)
            assert name in META['local'],\
            IO( f"""\nMaterial {name} not in the local bank {META['local']}\n
                Please refer to Documentation Material section
                https://pymiesim.readthedocs.io/en/latest/Material.html""")


        self.LocalDir =  META['local'][name]


    def _Plot(self):
        fig = plt.figure(figsize=(6,3.5))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'Wavelength $\lambda$ [m]')

        ax.plot(self.Data['wl0'], self.Data['n'], 'C0', label='n')
        ax.set_ylabel(r'Refractive index n')

        if 'wl1' in self.Data:
            ax1 = ax.twinx()
            ax1.plot(self.Data['wl1'], self.Data['k'], 'C1', label='k')
            ax1.set_ylabel(r'Extinction factor $\kappa$')

        ax.grid()
        ax.legend()
        ax1.legend()
        plt.tight_layout()



    def Plot(self):
        self._Plot()
        plt.show()

    def SaveFig(self, dir):
        dir = os.path.join(ZeroPath, 'docs/images', dir) + '.png'
        self._Plot()
        plt.savefig(dir)
        plt.close('all')


    def __repr__(self):
        return self.__name__

    def __str__(self):
        return self.__name__



















# --
