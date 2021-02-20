#!/usr/bin/env python

# -*- coding: utf-8 -*-

from unittest import TestCase
from numpy import linspace

from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.Source import PlaneWave
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.Sets import ScattererSet, ExperimentalSet

LightSource = PlaneWave(Wavelength = 450e-9, Polarization = 0)
Scat        = Sphere(Diameter = 300e-9, Index = 1.4, Source = LightSource)
Samp        = WMSample(g = 0.8, lc = 1e-5, D = 2.5, Nc = 1e4, Source = LightSource)
Detector    = LPmode(Mode = (0, 1,'h'), Sampling = 11, NA = 0.2)
Detector1   = Photodiode(Sampling = 11, NA = 0.2)
ScatSet     = ScattererSet(DiameterList = linspace(100e-9, 4500e-9, 11), RIList = 1.5, Source = LightSource)

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

        print('Test 0: passed')

    def test1(self):
        Detector = LPmode(Mode         = (1, 1,'h'),
                          Sampling     = 11,
                          NA           = 0.2,
                          GammaOffset  = 0,
                          PhiOffset    = 0,
                          CouplingMode = 'Centered')
        print('Test 1: passed')

    def test2(self):
        Scat.S1S2(Num=10)
        print('Test 2: passed')

    def test3(self):
        Scat.Field(Num=10)
        print('Test 3: passed')

    def test4(self):
        Scat.SPF(Num=10)
        print('Test 4: passed')

    def test5(self):
        Detector.Coupling(Scatterer = Scat)
        Detector1.Coupling(Scatterer = Scat)
        print('Test 5: passed')

    def test6(self):
        ScatSet = ScattererSet(DiameterList  = linspace(100e-9, 4500e-9, 11),
                               RIList        = 1.5,
                               Source        = LightSource)

        ScatSet.Qsca()
        print('Test 6: passed')

    def test7(self):
        Set = ExperimentalSet(ScattererSet = ScatSet,  Detectors = Detector)
        Set.DataFrame
        Set.Coupling
        print('Test 7: passed')

    def test8(self):
        Sample = WMSample(g      = 0.8,
                          lc     = 4e-5,
                          D      = 3/2,
                          Nc     = 1e4,
                          Source = LightSource)
        print('Test 8: passed')

    def test9(self):
        ScatSet = SampleSet(gList    = [0.8, 0.9],
                            LcList   = [1e-5, 2e-5],
                            D        = 3/2,
                            Nc       = 1e4,
                            Detector = Detector,
                            Source   = LightSource,
                            Npts     = 201)


        print('Test 9: passed')

if __name__ == '__main__':
    test = PrintingTest()
    test.Run()
