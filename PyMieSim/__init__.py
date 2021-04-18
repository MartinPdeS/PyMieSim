import json
import os
import csv
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

PATH = os.path.join( Path(__file__).parent, '_Material' )

from PyMieSim.BaseClasses import BaseMaterial

class Material(BaseMaterial):
    def __init__(self, name):
        with open(os.path.join(PATH, 'Meta.json'), 'r+' ) as f:
            META                     = json.load(f)
            assert name in META['local']

        self.dir =  META['local'][name]

    def Plot(self):
        data = self.LoadLocal(self.dir)
        print(data)
        plt.figure()
        plt.plot(data)
