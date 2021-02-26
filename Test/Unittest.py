#!/usr/bin/env python
# -*- coding: utf-8 -*-

from unittest import TestCase
from numpy import linspace, pi

from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.Source import PlaneWave, GaussianBeam
from PyMieSim.GLMT.python.GLMT import SPF
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.Sets import ScattererSet, ExperimentalSet, SampleSet

LightSource = PlaneWave(Wavelength = 450e-9, Polarization = 0)
Scat        = Sphere(Diameter = 300e-9, Index = 1.4, Source = LightSource)
Samp        = WMSample(g = 0.8, lc = 1e-5, D = 2.5, Nc = 1e4, Source = LightSource)
Detector    = LPmode(Mode = (0, 1,'h'), Sampling = 11, NA = 0.2)
Detector1   = Photodiode(Sampling = 11, NA = 0.2)
ScatSet     = ScattererSet(DiameterList = linspace(100e-9, 4500e-9, 11), RIList = 1.5, Source = LightSource)
phi         = linspace(-pi/2, pi/2,4)
theta       = linspace(-pi, pi,4)

class PrintingTest(TestCase):

    def Run(self):
        self.test00()
        self.test01()
        self.test02()
        self.test03()
        self.test04()
        self.test05()
        self.test06()
        self.test07()
        self.test08()
        self.test09()
        self.test10()

    def test00(self):
        Detector = LPmode(Mode         = (1, 1,'h'),
                          Sampling     = 11,
                          NA           = 0.2,
                          GammaOffset  = 0,
                          PhiOffset    = 0,
                          CouplingMode = 'Centered')

        print('Test 0: passed')

    def test01(self):
        Detector = LPmode(Mode         = (1, 1,'h'),
                          Sampling     = 11,
                          NA           = 0.2,
                          GammaOffset  = 0,
                          PhiOffset    = 0,
                          CouplingMode = 'Centered')
        print('Test 1: passed')

    def test02(self):
        Scat.S1S2(Num=10)
        print('Test 2: passed')

    def test03(self):
        Scat.FarField(Num=10)
        print('Test 3: passed')

    def test04(self):
        Scat.SPF(Num=10)
        print('Test 4: passed')

    def test05(self):
        Detector.Coupling(Scatterer = Scat)
        Detector1.Coupling(Scatterer = Scat)

        print('Test 5: passed')

    def test06(self):
        ScatSet = ScattererSet(DiameterList  = linspace(100e-9, 4500e-9, 11),
                               RIList        = 1.5,
                               Source        = LightSource)

        ScatSet.Qsca()

        print('Test 6: passed')

    def test07(self):
        Set = ExperimentalSet(ScattererSet = ScatSet,  Detectors = Detector)
        Set.DataFrame
        Set.Coupling

        print('Test 7: passed')

    def test08(self):
        Sample = WMSample(g      = 0.8,
                          lc     = 4e-5,
                          D      = 3/2,
                          Nc     = 1e4,
                          Source = LightSource)

        print('Test 8: passed')

    def test09(self):
        ScatSet = SampleSet(gList    = [0.8, 0.9],
                            LcList   = [1e-5, 2e-5],
                            D        = 3/2,
                            Nc       = 1e4,
                            Detector = Detector,
                            Source   = LightSource,
                            Npts     = 201)

        print('Test 9: passed')


    def test10(self):
        spf = SPF(Scat, LightSource, phi, theta)

        print('Test 10: passed')

    def test11(self):
        beam = GaussianBeam(Wavelength   = 0.632e-6,
                            NA           = 0.15,
                            Polarization = 0,
                            offset       = [5e-6,5e-6,5-6])

        beam.Amn(1,0)
        beam.Bmn(1,1)

if __name__ == '__main__':
    test = PrintingTest()
    test.Run()
























# -
