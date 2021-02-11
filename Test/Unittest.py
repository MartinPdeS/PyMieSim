#!/usr/bin/env python

# -*- coding: utf-8 -*-

from unittest import TestCase
from numpy import linspace

from PyMieSim.classes.Scattering import Scatterer
from PyMieSim.Physics import Source
from PyMieSim.classes.Detector import LPmode, Photodiode
from PyMieSim.classes.Sets import ScattererSet, ExperimentalSet

LightSource = Source(Wavelength = 450e-9, Polarization = 0, Power = 1,  Radius = 1)
Scat = Scatterer(Diameter = 300e-9, Source = LightSource, Index = 1.4)
Detector = LPmode(Mode = (0, 1,'h'), Sampling = 11, NA = 0.2)
ScatSet = ScattererSet(DiameterList = linspace(100e-9, 4500e-9, 11), RIList = 1.5, Source = LightSource)

class PrintingTest(TestCase):

    def test0(self):
        Detector = LPmode(Mode         = (1, 1,'h'),
                          Sampling     = 11,
                          NA           = 0.2,
                          GammaOffset  = 0,
                          PhiOffset    = 0,
                          CouplingMode = 'Centered')

    def test1(self):
        Detector = LPmode(Mode         = (1, 1,'h'),
                          Sampling     = 11,
                          NA           = 0.2,
                          GammaOffset  = 0,
                          PhiOffset    = 0,
                          CouplingMode = 'Centered')

    def test2(self):
         Scat.S1S2(Num=10)

    def test3(self):
        Scat.Field(Num=10)

    def test4(self):
        Scat.SPF(Num=10)

    def test5(self):
        Detector.Coupling(Scatterer = Scat)

    def test6(self):
        ScatSet = ScattererSet(DiameterList  = linspace(100e-9, 4500e-9, 11),
                               RIList        = 1.5,
                               Source        = LightSource)

        ScatSet.Qsca()

    def test7(self):
        Set = ExperimentalSet(ScattererSet = ScatSet,  Detectors = Detector)

        Set.DataFrame
