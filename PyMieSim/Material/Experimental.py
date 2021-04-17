import numpy as np
import os
from PyMieSim.Material.utils import CSV2Numpy, LoadCsvLocal, NearestIndex, LoadCsvOnline
from pathlib import Path


URL = { 'BK7'            : 'https://refractiveindex.info/data_csv.php?datafile=data/glass/schott/N-BK7.yml',
        'FusedSilica'    : 'https://refractiveindex.info/data_csv.php?datafile=data/main/SiO2/Malitson.yml',
         'SodaLimeGlass' : 'https://refractiveindex.info/data_csv.php?datafile=data/glass/misc/soda-lime/Rubin-clear.yml'
        }

PATH = Path(__file__).parent.parent

LOCAL = { 'BK7'           : os.path.join(PATH, '_Material/N-BK7.csv'),
          'FusedSilica'   : os.path.join(PATH, '_Material/FusedSilica.csv'),
          'SodaLimeGlass' : os.path.join(PATH, '_Material/SodaLimeGlass.csv'),
          'Silver'        : os.path.join(PATH, '_Material/Silver.csv'),
        }


def BK7Glass(wavelength):
    #CSV = LoadCsvOnline(URL['BK7'])
    CSV = LoadCsvLocal(LOCAL['BK7'])

    Boundary = [CSV[0][0], CSV[0][-1]]

    idx = NearestIndex( CSV[0], wavelength * 1e3 )

    return CSV[1][idx]


def FusedSilica(wavelength):
    #CSV = LoadCsvOnline(URL['FusedSilica'])
    CSV = LoadCsvLocal(LOCAL['FusedSilica'])

    Boundary = [CSV[0][0], CSV[0][-1]]

    idx = NearestIndex( CSV[0], wavelength * 1e3 )

    return CSV[1][idx]



def SodaLimeGlass(wavelength):
    #CSV = LoadCsvOnline(URL['SodaLimeGlass'])
    CSV = LoadCsvLocal(LOCAL['SodaLimeGlass'])

    Boundary = [CSV[0][0], CSV[0][-1]]

    idx = NearestIndex( CSV[0], wavelength * 1e3 )

    return CSV[1][idx]



#-
