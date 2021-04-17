#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from mayavi import mlab
from unittest import TestCase
from numpy import linspace, pi
import numpy as np
import unittest

from PyMieSim.Scatterer import Sphere, Cylinder, WMSample
from PyMieSim.Source import PlaneWave, GaussianBeam
from PyMieSim.GLMT.python.Sphere import SPF
from PyMieSim.Detector import LPmode, Photodiode, _Photodiode
from PyMieSim.Experiment import ScatSet, Setup, SourceSet, SampleSet
from PyMieSim.Mesh import FibonacciMesh
from PyMieSim.Plots import *
from unittest.mock import patch
from PyMieSim.Representations import S1S2


LightSource = PlaneWave(Wavelength = 450e-9, Polarization = 0)
Scat        = Sphere(Diameter = 300e-9, Index = 1.4, Source = LightSource)
Samp        = WMSample(g = 0.8, lc = 1e-5, D = 2.5, Nc = 1e4, Source = LightSource)
Detector    = LPmode(Mode = (0, 1,'h'), Sampling = 11, NA = 0.2)
Detector1   = Photodiode(Sampling = 11, NA = 0.2)
phi         = linspace(-pi/2, pi/2,4)
theta       = linspace(-pi, pi,4)





class DetectorTestCase(unittest.TestCase):

    def test00(self):
        from PyMieSim.Detector import Photodiode
        global Photodiode
        Photodiode = Photodiode(Sampling     = 11,
                                     NA           = 0.2,
                                     GammaOffset  = 0,
                                     PhiOffset    = 0,
                                     CouplingMode = 'Centered',
                                     Testing      = True)

        print('Photodiode initialisation passed')


    def test01(self):
        from PyMieSim.Detector import LPmode
        global LPmode
        LPmode = LPmode(Mode         = (0,1),
                             Sampling     = 11,
                             NA           = 0.2,
                             GammaOffset  = 0,
                             PhiOffset    = 0,
                             CouplingMode = 'Centered',
                             Testing      = True)

        print('LPmode initialisation passed')


    def test02(self):
        from PyMieSim.Detector import IntegratingSphere
        global IntegratingSphere
        IntegratingSphere = IntegratingSphere(Testing = True)

        print('IntegratingSphere initialisation passed')


    def test03(self):
        Photodiode.Plot()
        print('Photodiode Plotting passed')


    def test04(self):
        LPmode.Plot()
        print('LPmode Plotting passed')


    def test05(self):
        IntegratingSphere.Plot()
        print('IntegratingSphere Plotting passed')


class ScattererTestCase(unittest.TestCase):

    def test00(self):
        from PyMieSim.Scatterer import Sphere
        global sScat
        sScat = Sphere(Diameter = 300e-9,
                            Index    = 1.4,
                            Source   = LightSource,
                            Testing  = True)

        print('Sphere initialisation passed')


    def test01(self):
        from PyMieSim.Scatterer import Cylinder
        global cScat
        cScat = Cylinder(Diameter = 300e-9,
                              Index    = 1.4,
                              Source   = LightSource,
                              Testing  = True)

        print('Cylinder initialisation passed')


    def test02(self):
        sScat.S1S2(Num=10)
        print('Spherical Scatterer <S1S2> compute passed')


    def test03(self):
        sScat.Stokes(Num=10)
        print('Spherical Scatterer <Stokes> compute passed')


    def test04(self):
        sScat.FarField(Num=10)
        print('Spherical Scatterer <FarFields> compute passed')


    def test05(self):
        sScat.SPF(Num=10)
        print('Spherical Scatterer <SPF> compute passed')


    def test06(self):
        Photodiode.Coupling(Scatterer = sScat)
        print('<Photodiode> coupling passed')

    def test07(self):
        LPmode.Coupling(Scatterer = sScat)
        print('<LPmode> coupling passed')


    def test08(self):
        IntegratingSphere.Coupling(Scatterer = sScat)
        print('<IntegratingSphere> coupling passed')


    def test09(self):
        sScat.S1S2(10).Plot()
        print('<S1S2> coupling passed')


    def test10(self):
        sScat.Stokes(10).Plot()
        print('<Stokes> coupling passed')


    def test11(self):
        sScat.FarField(10).Plot()
        print('<FarField> coupling passed')


    def test12(self):
        sScat.SPF(10).Plot()
        print('<SPF> coupling passed')




class ExperiementTestCase(unittest.TestCase):

    def test00(self):
        global sScatSet
        sScatSet = ScatSet(DiameterList = linspace(100e-9, 4500e-9, 11),
                           IndexList      = 1.5)
        print('<ScattererSet> initialisation passed')









"""

class PrintingTest(TestCase):

    def Run(self):
        self.test09()
        self.test10()
        self.test11()
        self.test12()
        self.test13()
        self.test14()
        self.test15()
        self.test16()
        self.test17()
        self.test18()
        self.test19()
        self.test20()
        self.test21()
        self.test22()
        self.test23()
        self.test24()
        self.test25()
        self.test26()
        self.test27()
        self.test28()
        self.test29()


    def test09(self):
        self.sScatSet = ScatSet(DiameterList = linspace(100e-9, 4500e-9, 11), IndexList = 1.5)
        print('Test 09: <ScattererSet> initialisation passed')


    def test10(self):
        self.source = SourceSet(WavelengthList = 400e-9, PolarizationList = 0, SourceType = 'PlaneWave')
        self.ExpSet = Setup(ScattererSet = self.sScatSet, SourceSet = self.source, DetectorSet = self.Photodiode)
        print('Test 10: <Setup> initialisation passed')


    def test11(self):
        self.ExpSet.Efficiencies()
        print('Test 11: <ScattererSet> Qsca passed')


    def test12(self):
        self.ExpSet.Coupling
        print('Test 12: <Setup> dataframe compute passed')


    def test13(self):
        self.WMSample = WMSample(g      = 0.8,
                                 lc     = 4e-5,
                                 D      = 3/2,
                                 Nc     = 1e4,
                                 Source = LightSource)

        print('Test 13: WM sample initialisation passed')


    def test14(self):
        self.SampleSet = SampleSet(gList    = [0.8, 0.9],
                                   LcList   = [1e-5, 2e-5],
                                   D        = 3/2,
                                   Nc       = 1e4,
                                   Detector = Detector,
                                   Source   = LightSource,
                                   Npts     = 201)

        print('Test 14: SampleSet initialisation passed')


    def test15(self):
        self.Gbeam = GaussianBeam(Wavelength   = 1.3e-6,
                            NA           = 0.6,
                            Polarization = 0,
                            Offset       = [0e-6,0e-6,0e-6])

        print('Test 15: GaussianBeam beam initialisation passed')


    def test16(self):
        self.Gbeam.GetBSC(MaxOrder=3, save=False, Sampling=100)
        print('Test 16: GaussianBeam beam BSC compute passed')


    def test17(self):
        self.PWbeam = PlaneWave(Wavelength = 0.632e-6, Polarization = 0)
        print('Test 17: PlaneWave beam initialisation passed')


    def test18(self):
        self.PWbeam.GetBSC(MaxOrder=1, save=False)
        print('Test 18: GaussianBeam beam BSC compute passed')


    def test19(self):
        self.Photodiode.Footprint(Scatterer=self.sScat, Num=10)
        print('Test 19: Photodiode footprint compute passed')


    def test20(self):
        self.LPMode.Footprint(Scatterer=self.sScat, Num=10)
        print('Test 20: LPmode footprint compute passed')


    def test21(self):
        Scat = Cylinder(Diameter = 300e-9, Index = 1.4, Source = self.Gbeam)
        Scat.SPF(Num=10)
        print('Test 21: SPF compute with GLMT')


    def test22(self):
        Mesh = FibonacciMesh(MaxAngle    = pi,
                             Sampling    = 1000,
                             PhiOffset   = 0,
                             GammaOffset = 0)

        val0 = Scat.CrossSection(Mesh)
        val1 = Scat.Qsca * Scat.Area
        Rerror = np.abs(val0-val1)/val0
        assert Rerror < 1e-2
        print('Test 22: Validation QSca - CrossSection passed')


    def test23(self):
        Detector1 = _Photodiode(Sampling = 500, NA = 2.0)
        val0      = Scat.EnergyFlow(Detector1.Mesh)
        val1      = Detector1.Coupling(Scat)
        error     = np.abs(val0-val1)/val0
        assert error < 1e-2
        print('Test 23: Validation EnergyFlow - coupling passed')


    def test24(self):
        self.sScat.Stokes(Num=10)
        print("Test 24: Scatterer <Stokes> compute passed")


    def test25(self):
        self.ExpSet.Coupling(AsType='numpy')
        print("Test 25: <Experiment> 'numpy' output passed")


    def test26(self):
        self.ExpSet.Coupling(AsType='pymiesim')
        print("Test 26: <Experiment> 'ndarray' output passed")


    def test27(self):
        self.ExpSet.Coupling(AsType='dataframe')
        print("Test 27: <Experiment> 'dataframe' output passed")


    def test28(self):
        self.ExpSet.Coupling(AsType='optimizer')
        print("Test 28: <Experiment> 'optimizer' output passed")


    def test29(self):
        with patch('matplotlib.pyplot') as p:
            a = self.sScat.S1S2(Num=10)
            assert a.Plot() == None
            print('Test 29: Scatterer <S1S2> plot passed')
"""

def suite():
  suite = unittest.TestSuite()
  suite.addTest(DetectorTestCase('test00'))
  suite.addTest(DetectorTestCase('test01'))
  suite.addTest(DetectorTestCase('test02'))
  suite.addTest(DetectorTestCase('test03'))
  suite.addTest(DetectorTestCase('test04'))
  suite.addTest(DetectorTestCase('test05'))

  suite.addTest(ScattererTestCase('test00'))
  suite.addTest(ScattererTestCase('test01'))
  suite.addTest(ScattererTestCase('test02'))
  suite.addTest(ScattererTestCase('test03'))
  suite.addTest(ScattererTestCase('test04'))
  suite.addTest(ScattererTestCase('test05'))
  suite.addTest(ScattererTestCase('test06'))
  suite.addTest(ScattererTestCase('test07'))
  suite.addTest(ScattererTestCase('test08'))
  suite.addTest(ScattererTestCase('test09'))
  suite.addTest(ScattererTestCase('test10'))
  suite.addTest(ScattererTestCase('test11'))
  suite.addTest(ScattererTestCase('test12'))

  suite.addTest(ExperiementTestCase('test00'))

  return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner(failfast=True)
    runner.run(suite())
























# -
