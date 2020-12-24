import numpy as np
import pandas as pd
from typing import Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Optimizer import OptArray
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.BaseClasses import BaseScatterer
from PyMieCoupling.cpp.S1S2 import GetS1S2
from PyMieCoupling.classes.DataFrame import DataFrameCPU





class ScattererSet(object):

    def __init__(self,
                 DiameterList:    list,
                 RIList:          list,
                 Detector:        Union[LPmode, Photodiode],
                 Source:          Source,
                 Mode:            str = 'Centered',
                 Npts:            int = 201,
                 ):

        self.DiameterList, self.RIList = DiameterList, RIList

        self.Detector, self.Source, self.Mode = Detector, Source, Mode

        self.Coupling = np.empty( [len(self.RIList), len(self.DiameterList)] )

        maxPhi, minPhi = np.max(Detector.Meshes.Phi.Degree), np.min(Detector.Meshes.Phi.Degree)

        self.PhiVector = np.linspace(-np.pi/2, +np.pi/2,201)

        self.PhiVectorDetector = np.linspace(minPhi, maxPhi, 101)


    def GetCouplingFrame(self, Filter: list = ['None'] ):

        if not isinstance(Filter, list): Filter = [Filter]

        MI = pd.MultiIndex.from_product([Filter, self.DiameterList, self.RIList],
                                        names=['Filter','Diameter','RI',])

        df = DataFrameCPU(index = MI, columns = ['Coupling'])

        df.attrs['Filter'] = Filter

        for nr, RI in enumerate( self.RIList ):

            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter    = Diameter,
                                 Index       = RI,
                                 Source      = self.Source,
                                 Meshes      = self.Detector.Meshes)

                for Polar in df.attrs['Filter']:
                    self.Detector.Filter = Polar

                    Coupling = self.Detector.Coupling(Scatterer    = Scat,
                                                      Mode         = self.Mode)

                    df.at[(Polar, Diameter, RI),'Coupling'] = Coupling


        df.Coupling = df.Coupling.astype(float)

        df['Mean'] = df.groupby(['Filter','Diameter']).Coupling.transform('mean')

        df['STD'] = df.groupby(['Filter','Diameter']).Coupling.transform('std')

        df.DetectorNane = self.Detector._name

        return df



    def GetS1S2(self, num=201, n=0):
        List = np.empty( [len(self.RIList), len(self.DiameterList), self.Detector.Npts], dtype=np.complex  )
        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):
                SizeParam =  2 * np.pi * Diameter/self.Source.Wavelength

                S1S2 = GetS1S2(RI, SizeParam, self.PhiVector);

                List[nr, nd,:] = S1S2[n+1]




    def GetCouplingArray(self, Filter='None'):

        temp = np.empty( [ len(self.RIList), len(self.DiameterList) ] )


        for nr, RI in enumerate(self.RIList):
            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter  = Diameter,
                                 Index     = RI,
                                 Source    = self.Source,
                                 Meshes    = self.Detector.Meshes
                                 )

                Coupling = self.Detector.Coupling(Scatterer    = Scat,
                                                  Mode         = self.Mode)

                temp[nr, nd] = Coupling

        return OptArray(temp)



    def Plot(self, y: str = 'S1'):

        if 'S1' in  y: data = self.GetS1S2(1)

        elif "S2" in y: data = self.GetS1S2(2)

        data = np.abs(data)**2

        if str(y) == "S1": self.Plot_S1(Meshes, data)

        if str(y) == 'STD::S1': self.Plot_STDS1(Meshes, data)

        if str(y) == 'S2': self.Plot_S2(Meshes, data)

        if str(y) == 'STD::S2': self.Plot_STDS2(Meshes, data)



    def Plot_STDS1(self, Meshes, y):

        STDDiameter = y.std(axis=1)

        STDRI = y.std(axis=0)

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        ax.grid()

        for nr, RI in enumerate(self.RIList):

            ax.plot(self.PhiVector,
                    STDDiameter[nr],
                    '--',
                    label="RI:{0:.2f}".format(RI)
                    )


        for nd, Diameter in enumerate(self.DiameterList):

            ax.plot(self.PhiVector,
                    STDRI[nd],
                    label="Diam.:{0:.2e}".format(Diameter)
                    )


        ax.fill_between(self.PhiVector,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.Meshes.Phi.Boundary.Degree[1]) ,
                        label = 'Detector',
                        color = 'green',
                        alpha = 0.5)

        ax.set_xlabel(r'Scattering Angle [degree]')

        ax.set_ylabel(r'Scattered light intensity: $STD(S1^2)$ ')

        ax.set_yscale('log')

        ax.tick_params(labelsize=8)

        plt.legend(fontsize=6)

        plt.show(block=False)



    def Plot_STDS2(self, Meshes, y):

        STDDiameter = y.std(axis=1)

        STDRI = y.std(axis=0)

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        ax.grid()

        for nr, RI in enumerate(self.RIList):

            ax.plot(self.PhiVector,
                    STDDiameter[nr],
                    '--',
                    label="RI:{0:.2f}".format(RI)
                    )


        for nd, Diameter in enumerate(self.DiameterList):

            ax.plot(self.PhiVector,
                    STDRI[nd],
                    label="Diam.:{0:.2e}".format(Diameter)
                    )


        ax.fill_between(self.PhiVectorDetector,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.Meshes.Phi.Boundary.Degree[1]) ,
                        label = 'Detector',
                        color = 'green',
                        alpha = 0.5)

        ax.set_xlabel(r'Scattering Angle [degree]')

        ax.set_ylabel(r'Scattered light intensity: $STD(S2^2)$ ')

        ax.set_yscale('log')

        ax.tick_params(labelsize=8)

        plt.legend(fontsize=6)

        plt.show(block=False)



    def Plot_S1(self, Meshes, y):

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):

                plt.plot(self.PhiVector,
                         y[nr, nd],
                         label="RI:{0:.2f}; Diam.: {1:.3e}".format(RI, Diameter))


        ax.fill_between(self.PhiVectorDetector,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.Meshes.Phi.Boundary.Degree[1]) ,
                        label = 'Detector',
                        color = 'green',
                        alpha = 0.5)

        ax.grid()

        ax.set_xlabel(r'Scattering Angle [degree]')

        ax.set_ylabel(r'Scattered light intensity: $S1^2$ ')

        ax.set_yscale('log')

        ax.tick_params(labelsize=8)

        plt.legend(fontsize=6)

        plt.show(block=False)



    def Plot_S2(self, y):

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):

                plt.plot(self.PhiVector,
                         y[nr, nd],
                         label="RI:{0:.2f}; Diam.: {1:.3e}".format(RI, Diameter))


        ax.fill_between(self.PhiVectorDetector,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.Meshes.Phi.Boundary.Degree[1]) ,
                        label = 'Detector',
                        color = 'green',
                        alpha = 0.5)

        ax.grid()

        ax.set_xlabel(r'Scattering Angle [degree]')

        ax.set_ylabel(r'Scattered light intensity: $S2^2$ ')

        ax.set_yscale('log')

        ax.tick_params(labelsize=8)

        plt.legend(fontsize=6)

        plt.show(block=False)




class Scatterer(BaseScatterer):

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Meshes:      AngleMeshes  = None,
                 Acceptance:  list         = 20,
                 Samples:     int          = 1000,
                 GammaOffset: float        = 0,
                 PhiOffset:   float        = 0) -> None:

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.Acceptance = np.deg2rad(Acceptance)

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        if Meshes:
            self.Meshes = Meshes
        else:
            self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                      Samples     = Samples,
                                      PhiOffset   = PhiOffset,
                                      GammaOffset = GammaOffset)
                     




class FullScatterer(BaseScatterer):

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Samples:     int     = 1000):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.Acceptance = np.deg2rad(180)

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                  Samples     = Samples,
                                  PhiOffset   = 0,
                                  GammaOffset = 0)








# -
