import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt

from PyMieSim.Plots import PlotStructuredAmplitude
from PyMieSim.Plots import PlotStructuredAbs
from PyMieSim.utils import Direct2spherical, AngleUnit2DirectUnit



class Stokes(dict):
    def __array_finalize__(self, obj):
        pass

    def __init__(self, Field):
        pass



class SPF(dict):

    def __init__(self, Parent, Num=100, Distance=1.):

        Phi, Theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)
        self['Phi'], self['Theta'] = np.meshgrid(Phi, Theta)

        self['EPhi'], self['ETheta'] = Parent.Bind.sFields(Phi = Phi, Theta=Theta, R=1.)

        self['SPF'] = np.sqrt( self['EPhi'].__abs__()**2 + self['ETheta'].__abs__()**2 )


    def Plot(self, hold=False):
        Name = 'Scattering phase function'

        PlotStructuredAbs(self['SPF'],
                          self['Phi'],
                          self['Theta'],
                          Name='Scattering phase function')

        if hold: return
        else: mlab.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            SPF, EPhi, ETheta, Theta, Phi
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""



class S1S2(dict):

    def __init__(self, Parent, Num):

        self['Phi'] = np.linspace(0, 2*np.pi, Num)

        Phi = np.linspace(-np.pi,np.pi,Num);

        self['S1'], self['S2'] = Parent.Bind.S1S2(Phi = Phi)


    def Plot(self, hold=False):

        S1 = np.abs(self['S1'])
        S2 = np.abs(self['S2'])

        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 subplot_kw = {'projection':'polar'})

        axes[0].set_title('S1 function'); axes[1].set_title('S2 function')

        axes[0].plot(self['Phi'],
                     S1,
                     color = 'k')

        axes[0].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S1,
                             color = 'C0',
                             alpha = 0.4)


        axes[1].plot(self['Phi'],
                     S2,
                     color = 'k')

        axes[1].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S2,
                             color = 'C1',
                             alpha = 0.4)

        if hold: return
        else: plt.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            S1, S2, Phi
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""




class ScalarFarField(dict):

    def __init__(self, Num = 200, Parent = None, Distance=1.):

        Phi, Theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)

        self['Phi'], self['Theta'] = np.meshgrid(Phi, Theta)

        self['EPhi'], self['ETheta'] = Parent.Bind.sFields(Phi = Phi, Theta=Theta, R=1.)

        self['Distance'] = Distance

        self['SPF'] = np.sqrt( self['EPhi'].__abs__()**2 + self['ETheta'].__abs__()**2 )


    def Plot(self, hold=False):
        """Method plots the scattered Far-Field
        :math:`E_{\\phi}(\\phi,\\theta)^2 , E_{\\theta}(\\phi,\\theta)^2`.

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        """

        PlotStructuredAmplitude(self['EPhi'],
                                self['Phi'],
                                self['Theta'],
                                Name = 'E phi')

        PlotStructuredAmplitude(self['ETheta'],
                                self['Phi'],
                                self['Theta'],
                                Name = 'E theta')

        if hold: return
        else: mlab.show()




    def __repr__(self):

        return f"""
        Object:          Dictionary
        Keys:            EPhi, ETheta, Theta, Phi, Distance
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Theta'].shape}"""





class Footprint(dict):

    def __init__(self, Scatterer, Detector, Num=100):

        x, y = np.mgrid[-50: 50: complex(Num), -50: 50: complex(Num)]

        MaxAngle = np.abs( np.pi/2 - Detector.Mesh.Phi.Radian.min() )

        Phi, Theta = Direct2spherical(X=x, Y=y, MaxAngle=MaxAngle)

        Direct = AngleUnit2DirectUnit(Phi, Scatterer.Source.k)

        FarFieldPara, FarFieldPerp = Scatterer.uS1S2(Phi.flatten(), Theta.flatten())

        Perp =  \
        Detector.StructuredFarField(Num=Num, SFactor=16) * FarFieldPerp.reshape(Theta.shape)

        Para = \
        Detector.StructuredFarField(Num=Num, SFactor=16) * FarFieldPara.reshape(Theta.shape)

        n = 5

        FourierPara = np.fft.ifft2(Para, s=[512*n, 512*n])

        FourierPara = np.fft.fftshift(FourierPara).__abs__()**2

        FourierPerp = np.fft.ifft2(Perp,  s=[512*n, 512*n])

        FourierPerp = np.fft.fftshift(FourierPerp).__abs__()**2

        Direct = np.linspace(np.min(Direct), np.max(Direct), FourierPerp.shape[0])

        self['Map'] = FourierPara + FourierPerp

        self['DirectX'] = Direct; self['DirectY'] = Direct



    def Plot(self, hold=False):

        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.contourf(self['DirectX']*1e6, self['DirectY']*1e6, self['Map'], cmap='gray')

        ax.set_xlabel(r'Offset distance in X-axis [$\mu$m]')

        ax.set_ylabel(r'Offset distance in Y-axis [$\mu$m]')

        ax.set_title('Scatterer Footprint')

        if hold: return
        else: plt.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            Map, DirectX, DirectY
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""



# -
