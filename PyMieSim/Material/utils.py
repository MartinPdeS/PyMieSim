import csv
import urllib.request
import codecs
import numpy as np
import pandas as pd


def LoadCsvOnline(url):
    ftpstream = urllib.request.urlopen(url)
    Array = pd.read_csv(ftpstream, delimiter=',').T.to_numpy()
    bound = np.where(Array == 'wl')
    return Array[:, bound[0][0]:bound[1][0]].astype(float)


def LoadCsvLocal(directory):
    with open(directory) as csv_file:
        CSV = pd.read_csv(csv_file, delimiter=',')
        return CSV.T.to_numpy()


def CSV2Numpy(CSV):
    return np.asarray(list(CSV)[1:]).astype('float')


def NearestIndex(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx
