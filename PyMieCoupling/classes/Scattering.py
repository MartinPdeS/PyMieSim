
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib
import matplotlib
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from typing import Tuple, Union


from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Optimizer import OptArray
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint
from PyMieCoupling.classes.BaseClasses import BaseScatterer



class DataFrameCPU(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)
        self.Filter = None
        self.ax = None


    @property
    def Parallel(self):
        return self.xs('Parallel')


    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')


    def Plot(self, y, Scale='Linear'):

        for Polar in self.attrs['Filter']:
            self._plot(y, Polar, Scale)


    def _plot(self, y, Filter, Scale):

        self.ax = self.xs(Filter).unstack(1).plot(y       = y,
                                                  grid    = True,
                                                  figsize = (8,3.5),
                                                  title   = r'[{0}] Filter: {1} [Degree]'.format(self.DetectorNane, Filter),
                                                  ylabel  = y,
                                                  xlabel  = r'Scatterer diameter [m]')

        self.ax.tick_params(labelsize='small')
        self.ax.legend(bbox_to_anchor=(1, 1), ncol=1)

        if Scale == 'Logarithmic':
            self.ax.set_yscale('log')

        plt.subplots_adjust(right=0.8,)

        plt.show(block=False)




class ScattererSet(object):

    def __init__(self,
                 DiameterList:    list,
                 RIList:          list,
                 Detector:        Union[LPmode, Photodiode],
                 Source:          Source,
                 Mode:            str = 'Centered'
                 ):

        self.DiameterList = DiameterList

        self.Detector = Detector

        self.RIList = RIList

        self.Source = Source

        self.Mode = Mode

        self.Coupling = np.empty( [len(self.RIList), len(self.DiameterList)] )

        self.S1List = np.empty( [len(self.RIList), len(self.DiameterList), self.Detector.Npts], dtype=np.complex  )

        self.S2List = np.empty( [len(self.RIList), len(self.DiameterList), self.Detector.Npts], dtype=np.complex  )



    def GetCouplingFrame(self, Filter: list = ['None'] ):

        if not isinstance(Filter, list):
            Filter = [Filter]


        MI = pd.MultiIndex.from_product([Filter, self.DiameterList, self.RIList],
                                        names=['Filter','Diameter','RI',])

        df = DataFrameCPU(index = MI, columns = ['Coupling'])

        df.attrs['Filter'] = Filter

        for nr, RI in enumerate( self.RIList ):

            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter    = Diameter,
                                 Index       = RI,
                                 Source      = self.Source,
                                 Meshes      = self.Detector.FarField.Meshes)

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



    def GetS1(self, Boundary: str = 'Detector'):

        if Boundary == 'Detector':
            Meshes = self.Detector.Meshes

        elif Boundary == 'Full':
            Meshes = AngleMeshes(ThetaBound = np.array([-180, 180], copy=False),
                                 PhiBound   = np.array([-90,90], copy=False),
                                 ThetaNpts  = self.Detector.FarField.Meshes.ThetaNpts,
                                 PhiNpts    = self.Detector.FarField.Meshes.PhiNpts)

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter  = Diameter,
                                 Index     = RI,
                                 Source    = self.Source,
                                 Meshes    = Meshes
                                 )

                self.S1List[nr, nd,:] = Scat.S1S2[0]


        return self.S1List, Meshes



    def GetS2(self, Boundary: str = 'Detector'):

        if Boundary == 'Detector':
            Meshes = self.Detector.Meshes

        elif Boundary == 'Full':
            Meshes = AngleMeshes(ThetaBound = np.array([-180, 180], copy=False),
                                 PhiBound   = np.array([-90,90], copy=False),
                                 Npts       = self.Detector.Npts)

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter  = Diameter,
                                 Index     = RI,
                                 Source    = self.Source,
                                 Meshes    = Meshes
                                 )

                self.S2List[nr, nd,:] = Scat.S1S2[1]


        return self.S2List, Meshes



    def GetCouplingArray(self, Filter='None'):

        temp = np.empty( [ len(self.RIList), len(self.DiameterList) ] )


        for nr, RI in enumerate(self.RIList):
            for nd, Diameter in enumerate(self.DiameterList):

                Scat = Scatterer(Diameter  = Diameter,
                                 Index     = RI,
                                 Source    = self.Source,
                                 Meshes    = self.Detector.FarField.Meshes
                                 )

                Coupling = self.Detector.Coupling(Scatterer    = Scat,
                                                  Mode         = self.Mode)

                temp[nr, nd] = Coupling

        return OptArray(temp)



    def Plot(self, y: str = 'S1'):

        if 'S1' in  y:
            data, Meshes = self.GetS1('Full')

        elif "S2" in y:
            data, Meshes = self.GetS2('Full')

        data = np.abs(data)**2

        if str(y) == "S1":

            self.Plot_S1(Meshes, data)

        if str(y) == 'STD::S1':
            self.Plot_STDS1(Meshes, data)

        if str(y) == 'S2':
            self.Plot_S2(Meshes, data)

        if str(y) == 'STD::S2':
            self.Plot_STDS2(Meshes, data)



    def Plot_STDS1(self, Meshes, y):

        STDDiameter = y.std(axis=1)

        STDRI = y.std(axis=0)

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        ax.grid()

        for nr, RI in enumerate(self.RIList):

            ax.plot(Meshes.Phi.Vector.Degree,
                    STDDiameter[nr],
                    '--',
                    label="RI:{0:.2f}".format(RI)
                    )


        for nd, Diameter in enumerate(self.DiameterList):

            ax.plot(Meshes.Phi.Vector.Degree,
                    STDRI[nd],
                    label="Diam.:{0:.2e}".format(Diameter)
                    )


        ax.fill_between(self.Detector.FarField.Meshes.Phi.Vector.Degree,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.FarField.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.FarField.Meshes.Phi.Boundary.Degree[1]) ,
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

            ax.plot(Meshes.Phi.Vector.Degree,
                    STDDiameter[nr],
                    '--',
                    label="RI:{0:.2f}".format(RI)
                    )


        for nd, Diameter in enumerate(self.DiameterList):

            ax.plot(Meshes.Phi.Vector.Degree,
                    STDRI[nd],
                    label="Diam.:{0:.2e}".format(Diameter)
                    )


        ax.fill_between(self.Detector.FarField.Meshes.Phi.Vector.Degree,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.FarField.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.FarField.Meshes.Phi.Boundary.Degree[1]) ,
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

                plt.plot(Meshes.Phi.Vector.Degree,
                         y[nr, nd],
                         label="RI:{0:.2f}; Diam.: {1:.3e}".format(RI, Diameter))


        ax.fill_between(self.Detector.FarField.Meshes.Phi.Vector.Degree,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.FarField.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.FarField.Meshes.Phi.Boundary.Degree[1]) ,
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



    def Plot_S2(self, Meshes, y):

        fig = plt.figure(figsize=(7,3))

        ax = fig.add_subplot(1,1,1)

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):

                plt.plot(Meshes.Phi.Vector.Degree,
                         y[nr, nd],
                         label="RI:{0:.2f}; Diam.: {1:.3e}".format(RI, Diameter))


        ax.fill_between(self.Detector.FarField.Meshes.Phi.Vector.Degree,
                        ax.get_ylim()[0],
                        ax.get_ylim()[1]*3,
                        #where = (Meshes.Phi.Vector.Degree > self.Detector.FarField.Meshes.Phi.Boundary.Degree[0]) & (Meshes.Phi.Vector.Degree < self.Detector.FarField.Meshes.Phi.Boundary.Degree[1]) ,
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
    """Object containing all scatterer-related attributes.

    Parameters
    ----------
    diameter : float
        Diameter of the scatterer.
    wavelength : float
        Wavelength of the incident lightfield.
    index : float
        Refractive index of the scatterer.
    npts : int
        Number of points for the full solid angle of the far-field, later to
        be interpolated.

    Attributes
    ----------
    Full : <Fields class>
        It represents the entire Far-field representation of the scatterer.
    ComputeS1S2 : type
        Methode using package PyMieScatt to compute S1 and S2 parameter form mu value.
    diameter
    wavelength
    index
    npts

    """

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Npts:        int         = None,
                 Meshes:      AngleMeshes  = None,
                 ThetaBound:  list        = [-180, 180],
                 ThetaOffset: float       = 0,
                 PhiBound:    list        = [-180, 180],
                 PhiOffset:   float       = 0) -> None:

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._FarField, self._S1S2 = (None,)*4

        if Meshes:
            self.Meshes = Meshes
        else:
            self.Meshes = AngleMeshes(ThetaBound  = np.asarray(ThetaBound),
                                      PhiBound    = np.asarray(PhiBound),
                                      ThetaNpts   = Npts,
                                      PhiNpts     = Npts,
                                      PhiOffset   = PhiOffset,
                                      ThetaOffset = ThetaOffset)




class FullScatterer(BaseScatterer):
    """Object containing all scatterer-related attributes.

    Parameters
    ----------
    diameter : float
        Diameter of the scatterer.
    wavelength : float
        Wavelength of the incident lightfield.
    index : float
        Refractive index of the scatterer.
    npts : int
        Number of points for the full solid angle of the far-field, later to
        be interpolated.

    Attributes
    ----------
    Full : <Fields class>
        It represents the entire Far-field representation of the scatterer.
    ComputeS1S2 : type
        Methode using package PyMieScatt to compute S1 and S2 parameter form mu value.
    diameter
    wavelength
    index
    npts

    """

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Npts:        int    = 201):

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        self._Stokes, self._SPF, self._FarField, self._S1S2 = (None,)*4

        self.Meshes = AngleMeshes(ThetaBound  = np.asarray([-180,180]),
                                  PhiBound    = np.asarray([0,180]),
                                  ThetaNpts   = Npts,
                                  PhiNpts     = Npts,
                                  PhiOffset   = 0,
                                  ThetaOffset = 0)








# -
