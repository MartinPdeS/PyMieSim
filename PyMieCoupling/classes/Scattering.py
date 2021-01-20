import numpy as np
import pandas as pd
from typing import Tuple, Union
from scipy.special import gamma
import matplotlib

import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
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
                 Sampling:     int          = 1000,
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
                                      Sampling     = Sampling,
                                      PhiOffset   = PhiOffset,
                                      GammaOffset = GammaOffset)





class FullScatterer(BaseScatterer):

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Sampling:     int     = 1000):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.Acceptance = np.deg2rad(180)

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._Parallel, self._Perpendicular, self._S1S2 = (None,)*5

        self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                  Sampling     = Sampling,
                                  PhiOffset   = 0,
                                  GammaOffset = 0)



class SampleSet(object):

    def __init__(self,
                 gList:           list,
                 LcList:          list,
                 D:               float,
                 Nc:              float,
                 Detector:        Union[LPmode, Photodiode],
                 Source:          Source,
                 Mode:            str = 'Centered',
                 Npts:            int = 201,
                 ):

        self.gList, self.LcList = gList, LcList

        self.D = D; self.Nc = Nc

        self.Source = Source

        self.Detector, self.Source, self.Mode = Detector, Source, Mode

        self.Coupling = np.empty( [len(self.LcList), len(self.gList)] )

        maxPhi, minPhi = np.max(Detector.Meshes.Phi.Degree), np.min(Detector.Meshes.Phi.Degree)

        self.PhiVector = np.linspace(-np.pi/2, np.pi/2,201)

        self.PhiVectorDetector = np.linspace(minPhi, maxPhi, 101)


    def GetCouplingFrame(self, Filter: list = ['None'] ):

        if not isinstance(Filter, list): Filter = [Filter]

        MI = pd.MultiIndex.from_product([Filter, self.gList, self.LcList], names=['Filter','g','lc',])

        df = DataFrameCPU(index = MI, columns = ['Coupling'])

        df.attrs['Filter'] = Filter

        for ng, g in enumerate( self.gList ):

            for nlc, lc in enumerate(self.LcList):

                Scat = Sample(g           = g,
                              lc          = lc,
                              D           = self.D,
                              Nc          = self.Nc,
                              Source      = LightSource,
                              Meshes      = self.Detector.Meshes)

                for Polar in df.attrs['Filter']:
                    self.Detector.Filter = Polar

                    Coupling = self.Detector.Coupling(Scatterer = Scat, Mode = self.Mode)

                    df.at[(Polar, g, lc),'Coupling'] = Coupling


        df.Coupling = df.Coupling.astype(float)

        df.DetectorNane = self.Detector._name

        return df




    def GetCouplingArray(self, Filter='None'):

        temp = np.empty( [ len(self.gList), len(self.LcList) ] )


        for ng, g in enumerate( self.gList ):

            for nlc, lc in enumerate(self.LcList):

                Samp = Sample(g           = g,
                              lc          = lc,
                              D           = self.D,
                              Nc          = self.Nc,
                              Source      = self.Source,
                              Meshes      = self.Detector.Meshes)

                Coupling = self.Detector.Coupling(Scatterer = Samp, Mode = self.Mode)

                temp[ng, nlc] = Coupling

        return OptArray(temp)



    def Plot(self, y: str = 'S1'):

        pass



class Sample(object):
    def __init__(self,
                 g,
                 lc,
                 D,
                 Nc,
                 Source,
                 Meshes      = None,
                 Acceptance  = 20,
                 Sampling    = 1000,
                 GammaOffset = 0,
                 PhiOffset   = 0) -> None:

        self.g  = g
        self.lc = lc
        self.D  = D
        self.Nc = Nc
        self.Acceptance = Acceptance; self.Source = Source

        if Meshes:
            self.Meshes = Meshes
        else:
            self.Meshes = AngleMeshes(MaxAngle    = self.Acceptance,
                                      Sampling    = Sampling,
                                      PhiOffset   = PhiOffset,
                                      GammaOffset = GammaOffset)


        self.Parallel = self.GetScalar(self.Meshes.Theta.Radian, self.Meshes.Phi.Radian)

        self.Perpendicular = self.Parallel


    def GetScalar(self, Phi, Theta):

        k = self.Source.k

        term0 = 2 * self.Nc * self.lc * gamma(self.D/2) / np.sqrt(np.pi) * k**4

        term1 = (1-np.sin(Phi-np.pi/2)**2*np.cos(Theta + self.Source.Polarization.Radian)**2)

        term2 = (1 + (2* k * self.lc * np.sin((Phi-np.pi/2)/2)**2)**(self.D/2))

        return term0 * term1 / term2


    def Plot(self, num=200):
        import cartopy.crs as ccrs

        ThetaMesh, PhiMesh = np.mgrid[0:2*np.pi:complex(num), -np.pi/2:np.pi/2:complex(num)]

        scalar = self.GetScalar(PhiMesh, ThetaMesh+np.pi/2)

        n=0

        fig0, axes0 = plt.subplots(nrows      = 1,
                                   ncols      = 2,
                                   figsize    = (8, 4),
                                   subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea()}
                                  )

        im = axes0[n].contourf(
                                np.rad2deg(ThetaMesh),
                                np.rad2deg(PhiMesh),
                                scalar.real,
                                cmap = 'inferno',
                                shading='auto',
                                transform=ccrs.PlateCarree(),
                                #levels=50
                                )

        gl = axes0[n].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
        gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
        plt.colorbar(mappable=im, fraction=0.046, orientation='vertical', ax=axes0[n])
        axes0[n].set_title('Real Part Far-Field')
        axes0[n].set_ylabel(r'Angle $\phi$ [Degree]')
        axes0[n].set_xlabel(r'Angle $\theta$ [Degree]')


        n += 1
        im = axes0[n].contourf(
                                np.rad2deg(ThetaMesh),
                                np.rad2deg(PhiMesh),
                                scalar.imag,
                                cmap = 'inferno',
                                shading='auto',
                                transform=ccrs.PlateCarree(),
                                #levels=50
                                )

        gl = axes0[n].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
        gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
        plt.colorbar(mappable=im, fraction=0.046, orientation='vertical', ax=axes0[n])
        axes0[n].set_title('Imaginary Part Far-Field')
        axes0[n].set_ylabel(r'Angle $\phi$ [Degree]')
        axes0[n].set_xlabel(r'Angle $\theta$ [Degree]')



# -
