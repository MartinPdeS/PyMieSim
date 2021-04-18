#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import csv
import urllib.request
import pandas as pd
import numpy as np
from pathlib import Path

PATH = os.path.join( Path(__file__).parent )

from PyMieSim.BaseClasses import BaseMaterial



def LoadOnline(url):
    ftpstream = urllib.request.urlopen(url)
    Array = pd.read_csv(ftpstream, delimiter=',').T.to_numpy()
    bound = np.where(Array == 'wl')

    if len(bound) == 0:
        data = { 'wl0' : Array[0].astype(float),
                 'n'   : Array[1].astype(float)}

    if len(bound) == 1:
        bound = bound[0]
        data = { 'wl0' : Array[0].astype(float),
                 'n'   : Array[1].astype(float)}

    if len(bound) == 2:
        bound = (bound[0][0],bound[1][0])
        data = { 'wl0' : Array[0][:bound[1]].astype(float),
                 'n'   : Array[1][:bound[1]].astype(float),
                 'wl1' : Array[0][bound[1]+1:].astype(float),
                 'k'   : Array[1][bound[1]+1:].astype(float) }

    return data


def LoadOnlineSave(url, filename):
    dict_data = LoadOnline(url)
    directory = os.path.join( PATH, 'csv', filename + '.csv' )

    with open(directory, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['wl0', 'n', 'wl1', 'k'])
        writer.writeheader()
        writer.writerow(dict_data)

    with open(os.path.join(PATH, 'Meta.json'), 'r+' ) as f:
        META                     = json.load(f)
        META['remote'][filename] = url
        META['local'][filename]  = filename + '.csv'
        f.seek(0)
        json.dump(META, f, indent=4)

    print(META)

#-
