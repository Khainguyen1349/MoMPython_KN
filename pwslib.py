import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d.axes3d as Axes3D
#from matplotlib import cm, colors
from copy import deepcopy
from scipy import interpolate



class field:
    def __init__(self,E_field,X,Y,z,Npoints,k0,omega,eps_0,eps_c):
        self.DX = E_field[:,:,0]
        self.DY = E_field[:,:,1]
        self.DZ = E_field[:,:,2]
        self.X = X
        self.Y = Y
        self.Z = z
        self.N = Npoints
        self.D = np.swapaxes(np.stack([self.DX, self.DY, self.DZ],axis=1),1,2)

def spectral(E,Nspec,k0):
    dx = E.X[0,1] - E.X[0,0]
    dy = E.Y[1,0] - E.Y[0,0]
    kx = (-1/2 + (np.arange(Nspec))/Nspec)*(2*np.pi/dx)
    ky = (-1/2 + (np.arange(Nspec))/Nspec)*(2*np.pi/dy)
    KX,KY = np.meshgrid(kx,ky)
    KZ = np.conj(np.sqrt((k0**2-KX**2-KY**2)*(1+0j)))
        
    sE = deepcopy(E)
    sE.X = KX
    sE.Y = KY
    sE.Z = KZ
    sE.N = Nspec
    
    sE.DX = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(E.DX)))
    sE.DY = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(E.DY)))
    sE.DZ = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(E.DZ)))
    sE.D = np.swapaxes(np.stack([sE.DX, sE.DY, sE.DZ],axis=1),1,2)
    return sE

def propagation(E,Z,Nspec,k0):
    dx = E.X[0,1] - E.X[0,0]
    dy = E.Y[1,0] - E.Y[0,0]
    kx = (-1/2 + (np.arange(Nspec))/Nspec)*(2*np.pi/dx)
    ky = (-1/2 + (np.arange(Nspec))/Nspec)*(2*np.pi/dy)
    KX,KY = np.meshgrid(kx,ky)
    KZ = np.conj(np.sqrt((k0**2-KX**2-KY**2)*(1+0j)))

    propagator = np.fft.ifftshift(np.exp(-1j*(Z-E.Z)*KZ))
        
    E_propagated = deepcopy(E)
    E_propagated.Z = E.Z
        #Matlab : E_rec.D = ifft2(propagateur.*fft2(E_anal.D));
    E_propagated.DX = np.fft.ifft2(np.multiply(propagator,np.fft.fft2(E.DX)))
    E_propagated.DY = np.fft.ifft2(np.multiply(propagator,np.fft.fft2(E.DY)))
    E_propagated.DZ = np.fft.ifft2(np.multiply(propagator,np.fft.fft2(E.DZ)))
    E_propagated.D = np.swapaxes(np.stack([E_propagated.DX, E_propagated.DY, E_propagated.DZ],axis=1),1,2)
    return E_propagated

def displayE(E) :
    fig, axs = plt.subplots(3, 2)
    axs[0, 0].pcolormesh(E.X,E.Y,abs(E.DX),shading='auto')
    axs[0, 0].set_title('X-Y-Z component - magnitude')
    axs[0, 1].pcolormesh(E.X,E.Y,np.angle(E.DX),shading='auto')
    axs[0, 1].set_title('X-Y-Z component - phase')
    
    axs[1, 0].pcolormesh(E.X,E.Y,abs(E.DY),shading='auto')
    axs[1, 1].pcolormesh(E.X,E.Y,np.angle(E.DY),shading='auto')
    
    axs[2, 0].pcolormesh(E.X,E.Y,abs(E.DZ),shading='auto')
    axs[2, 1].pcolormesh(E.X,E.Y,np.angle(E.DZ),shading='auto')
    
    for ax in axs.flat:
        ax.set(xlabel='x-axis', ylabel='y-axis')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()