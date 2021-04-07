import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt

from PyMieSim.Plots import StructuredAmplitude, StokesPlot, StructuredAbs
from PyMieSim.utils import Direct2spherical, AngleUnit2DirectUnit


class Stokes(dict): # https://en.wikipedia.org/wiki/Stokes_parameters

    def __init__(self, Parent, Num=100, Distance=1.):

        self.Polarization = Parent.Source.Polarization.Radian

        Phi, Theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)

        self['Phi'], self['Theta'] = np.meshgrid(Phi, Theta)

        EPhi, ETheta = Parent.Bind.sFields(Phi = Phi, Theta=Theta, R=1.)

        self['I'] = np.abs(EPhi)**2 + np.abs(ETheta)**2

        self['Q'] = np.abs(EPhi)**2 - np.abs(ETheta)**2

        self['U'] = +2 * np.real(EPhi*ETheta.conjugate())

        self['V'] = -2 * np.imag(EPhi*ETheta.conjugate())


    def Plot(self):
        Name = 'Scattering phase function'

        StokesPlot(I            = self['I'],
                   Q            = self['Q'],
                   U            = self['U'],
                   V            = self['V'],
                   Phi          = self['Phi'],
                   Theta        = self['Theta'],
                   Name         = 'Stokes Parameter',
                   Polarization = self.Polarization)

    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            S1, S2, S3,, S4, Theta, Phi
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['S1'].shape}"""


class SPF(dict):

    def __init__(self, Parent, Num=100, Distance=1.):

        self.Polarization = Parent.Source.Polarization.Radian

        Phi, Theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)

        self['Phi'], self['Theta'] = np.meshgrid(Phi, Theta)

        self['EPhi'], self['ETheta'] = Parent.Bind.sFields(Phi = Phi, Theta=Theta, R=1.)

        self['SPF'] = np.sqrt( self['EPhi'].__abs__()**2 + self['ETheta'].__abs__()**2 )


    def Plot(self):
        Name = 'Scattering phase function'

        StructuredAbs(Scalar       = self['SPF'],
                      Phi          = self['Phi'],
                      Theta        = self['Theta'],
                      Name         = 'Scattering phase function',
                      Polarization = self.Polarization)



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


    def Plot(self):

        S1 = np.abs(self['S1'])
        S2 = np.abs(self['S2'])

        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 subplot_kw = {'projection':'polar'})

        axes[0].set_title('S1 function'); axes[1].set_title('S2 function')

        axes[0].plot(self['Phi'], S1,  color = 'k')

        axes[0].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S1,
                             color = 'C0',
                             alpha = 0.4)


        axes[1].plot(self['Phi'], S2, color = 'k')

        axes[1].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S2,
                             color = 'C1',
                             alpha = 0.4)


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            S1, S2, Phi
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""




class ScalarFarField(dict):

    def __init__(self, Num = 200, Parent = None, Distance=1.):

        self.Polarization = Parent.Source.Polarization.Radian

        Phi, Theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)

        self['Phi'], self['Theta'] = np.meshgrid(Phi, Theta)

        self['EPhi'], self['ETheta'] = Parent.Bind.sFields(Phi = Phi, Theta=Theta, R=1.)

        self['Distance'] = Distance

        self['SPF'] = np.sqrt( self['EPhi'].__abs__()**2 + self['ETheta'].__abs__()**2 )


    def Plot(self):
        """Method plots the scattered Far-Field
        :math:`E_{\\phi}(\\phi,\\theta)^2 , E_{\\theta}(\\phi,\\theta)^2`.

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        """

        StructuredAmplitude(Scalar       = self['EPhi'],
                            Phi          = self['Phi'],
                            Theta        = self['Theta'],
                            Name         = u'E_φ',
                            Polarization = self.Polarization)

        StructuredAmplitude(Scalar       = self['ETheta'],
                            Phi          = self['Phi'],
                            Theta        = self['Theta'],
                            Name         = u'E_θ',
                            Polarization = self.Polarization)



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
        Detector.FarField(Num=Num, Structured=True) * FarFieldPerp.reshape(Theta.shape)

        Para = \
        Detector.FarField(Num=Num, Structured=True) * FarFieldPara.reshape(Theta.shape)

        n = 5

        FourierPara = np.fft.ifft2(Para, s=[512*n, 512*n])

        FourierPara = np.fft.fftshift(FourierPara).__abs__()**2

        FourierPerp = np.fft.ifft2(Perp,  s=[512*n, 512*n])

        FourierPerp = np.fft.fftshift(FourierPerp).__abs__()**2

        Direct = np.linspace(np.min(Direct), np.max(Direct), FourierPerp.shape[0])

        self['Map'] = FourierPara + FourierPerp

        self['DirectX'] = Direct; self['DirectY'] = Direct



    def Plot(self):

        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.contourf(self['DirectX']*1e6, self['DirectY']*1e6, self['Map'], cmap='gray')

        ax.set_xlabel(r'Offset distance in X-axis [$\mu$m]')

        ax.set_ylabel(r'Offset distance in Y-axis [$\mu$m]')

        ax.set_title('Scatterer Footprint')



    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            Map, DirectX, DirectY
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""



# -
