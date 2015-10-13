import numpy as np
import pylab as pl
import scipy.fftpack as sfft
import scipy.ndimage.filters as sft

from genscript.extendclass import *
import genscript.myarray as ar
import genscript.fft as fft




def Poisson_3d(d, boxsize=None, return_hessian=False):
    s=d.shape
    ng = s[0]

    kmin = 2.*np.pi/float(ng)

    # ->> Do NOT shift freq for `sin(k)' kernel
    kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng/2+1].astype(float)
    k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)

    if boxsize!=None:
        k2*=(float(ng)/float(boxsize))**2.

    #kx = np.fromfunction(lambda x,y,z:x, (ng,ng,ng/2+1))
    #ky = np.fromfunction(lambda x,y,z:y, (ng,ng,ng/2+1))
    #kz = np.fromfunction(lambda x,y,z:z, (ng,ng,ng/2+1))
    #kx[np.where(kx > s[0]/2)] -= s[0]
    #ky[np.where(ky > s[0]/2)] -= s[0]
    #kz[np.where(kz > s[0]/2)] -= s[0]

    #k2=(kx**2.+ky**2.+kz**2.)*kmin**2.


    if return_hessian==False:

        dk=np.fft.rfftn(d)
        phik=-dk/k2
        phik[0,0,0]=0

        return np.fft.irfftn(phik)

    else:
        kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng].astype(float)
	ki = np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), np.sin(kz*kmin/2.)])*2.
        k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)
	
        if boxsize!=None:
            k2*=(float(ng)/float(boxsize))**2.
	    ki*=float(ng)/float(boxsize)

	print 'ki shape:', ki.shape

        dk=np.fft.fftn(d)
        phik=-dk/k2
        phik[0,0,0]=0

        phi = np.fft.ifftn(phik)
	phi_ij = np.zeros([len(s), len(s)]+list(d.shape)) 

        for i in range(3):
	    for j in range(3):
	        phi_ij[i,j] = np.fft.ifftn(-phik*ki[i]*ki[j])

	return phi.astype(float), phi_ij.astype(float)




def Poisson_1d(d, boxsize):
    s=d.shape
    ng = s[0]
    kmin = 2.*np.pi/float(ng)

    kx = np.arange(ng/2+1) #np.fromfunction(lambda x,y,z:x, (ng,ng,ng/2+1))
    kx[np.where(kx > s[0]/2)] -= s[0]

    k2=kx**2.*kmin**2.
    dk=np.fft.rfftn(d)

    phik=-dk/k2
    phik[0]=0

    return np.fft.irfftn(phik)




def Poisson(d, fft_axes=None):
    # ->> general Poisson solver for arbitrary dimension and data type.
    dk=np.fft.fftn(d, axes=fft_axes)

    s, ks= d.shape, dk.shape
    k=ar.mgrid(ks)

    #for i in range():


    return




def Laplacian(d, boxsize=None, Lap_type=1):

    if Lap_type==0:
        l=sft.laplace(d)
    elif Lap_type==1:
        l=-6*d
	for i in range(3):
            l+=np.roll(d, 1, axis=i)+np.roll(d, -1, axis=i)
    elif Lap_type==2:
        grd_phi=np.array([ np.gradient(np.gradient(d)[i]) for i in \
	         range(3) ])
	l=(grd_phi[0,0]+grd_phi[1,1]+grd_phi[2,2])

    if boxsize==None:
        dl=1.
    else:
        dl=boxsize/d.shape[0]

    return l/dl**2.




def Laplacian2d(d, boxsize=None):
    l=sft.laplace(d)

    if boxsize==None:
        dl=1.
    else:
        dl=boxsize/d.shape[0]

    return l/dl**2.



def Hessian(d, boxsize=None, numerical_type=0):

    dim=len(d.shape)
    print 'Hessian dimension:', dim

    if numerical_type==0:
        raise Exception


    if numerical_type==1:
        dij = np.zeros([dim, dim]+list(d.shape)) 

        for i in range(dim):
	    for j in range(dim):
	        dij[i,j] = np.roll(np.roll(d, 1, axis=i), 1, axis=j) - \
		           np.roll(d, 1, axis=i) - np.roll(d, 1, axis=j)\
			   + d
    elif numerical_type==2:
        dij=np.array([ np.gradient(np.gradient(d)[i]) for i in \
	               range(dim) ])
    else:
        raise Exception


    if boxsize==None:
        dl=1.
    else:
        dl=np.float(boxsize)/d.shape[0]

    return dij/dl**2.
