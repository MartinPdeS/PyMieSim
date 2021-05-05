#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import csv
import urllib.request
import pandas as pd
from socket import timeout
import numpy as np

from PyMieSim.BaseClasses import BaseMaterial
from PyMieSim.Directories import NPZPath, MaterialPath


def LoadOnline(url):
    try:
        ftpstream = urllib.request.urlopen(url)
    except ConnectionResetError:
        print("==> ConnectionResetError")
        pass
    except timeout:
        print("==> Timeout")
        pass
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


def LoadOnlineSave(url, filename, unit=1e-6):
    dict_data = LoadOnline(url)
    directory = os.path.join( NPZPath, filename )

    if 'wl1' not in dict_data:
         np.savez(directory,
                  wl0 = dict_data['wl0'] * unit,
                  n = dict_data['n'])
    else:
        np.savez(directory,
                 wl0 = dict_data['wl0'] * unit,
                 n = dict_data['n'],
                 wl1 = dict_data['wl1'] * unit,
                 k = dict_data['k']
                 )

    with open(os.path.join(MaterialPath, 'Meta.json'), 'r+' ) as f:
        META                     = json.load(f)
        META['remote'][filename] = url
        META['local'][filename]  = filename + '.npz'
        f.seek(0)
        json.dump(META, f, indent=4)


#-
