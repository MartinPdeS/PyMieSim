#!/usr/bin/env python

# -*- coding: utf-8 -*-

from unittest import TestCase
from numpy import linspace

from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.Physics import Source
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.Sets import ScattererSet, ExperimentalSet

LightSource = Source(Wavelength = 450e-9, Polarization = 0, Power = 1,  Radius = 1)
Scat = Sphere(Diameter = 300e-9, Source = LightSource, Index = 1.4)
Detector = LPmode(Mode = (0, 1,'h'), Sampling = 11, NA = 0.2)
ScatSet = ScattererSet(DiameterList = linspace(100e-9, 4500e-9, 11), RIList = 1.5, Source = LightSource)

class PrintingTest(TestCase):

    def Run(self):
        self.test0()
        self.test1()
        self.test2()
        self.test3()
        self.test4()
        self.test5()
        self.test6()
        self.test7()
        self.test8()

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

    def test8(self):
        Sample = WMSample(g      = 0.8,
                          lc     = 4e-5,
                          D      = 3/2,
                          Nc     = 1e4,
                          Source = LightSource)


if __name__ == '__main__':
    test = PrintingTest()
    test.Run()
