''' ->> modified version from 'poisson.py', originally from project `vgradth' <<- 
'''
import numpy as np
import pylab as pl
import scipy.fftpack as sfft
import scipy.ndimage.filters as sft

from genscript.extendclass import *
import genscript.myarray as ar
import genscript.fft as fft




def Poisson3d(d, boxsize=None, return_hessian=False, return_gradient=False, smooth_R=None, smooth_type=None):
    ''' ->> Poisson solver of given field 'd', also return gradient or hessian field <<- 
    '''

    ng = d.shape[0]
    ndim=len(d.shape)
    kmin = 2.*np.pi/float(ng)

    # ->> Fourier space arguments <<- #
    kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng].astype(float)

    ki = 2.*np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), np.sin(kz*kmin/2.)])
    k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)
    
    if boxsize!=None:
        k2*=(float(ng)/float(boxsize))**2.
        ki*=float(ng)/float(boxsize)
    
    print 'ki shape:', ki.shape


    # ->> if do_smoothing <<- #
    if (smooth_R!=None)&(smooth_type!=None):
        do_smooth=True

	if smooth_type=='Gaussian':
            W=np.exp(-k2*smooth_R**2./2.)
	else:
	    raise Exception
    else:
        do_smooth=False
	W=1.

    print 'do_smooth=', do_smooth
    
    # ->> Fourier transform & phi <<- #
    dk=np.fft.fftn(d)*W
    phik=-dk/k2
    phik[0,0,0]=0

    # ->> get Phi <<- #
    phi = np.fft.ifftn(phik)

    # ->> vector & tensors <<- #
    if return_hessian==True:
	phi_ij = np.zeros([ndim, ndim]+list(d.shape)) 

        for i in range(ndim):
	    for j in range(ndim):
	        phi_ij[i,j] = np.fft.ifftn(-phik*ki[i]*ki[j])

    if return_gradient==True:
	phi_i = np.zeros([ndim]+list(d.shape)) 

        for i in range(ndim):
            phi_i[i] = np.fft.ifftn(1.j*phik*ki[i])
	    #print 'phi_i real part:',np.min(phi_i[i].real) , np.max(phi_i[i].real)
	    #print 'phi_i imag part:',np.min(phi_i[i].imag) , np.max(phi_i[i].imag)

    # ->> return <<- #
    if (return_hessian==False)&(return_gradient==False):
        return phi.astype(float)

    elif (return_hessian==True)&(return_gradient==False):
        return phi.astype(float), phi_ij.astype(float)

    elif (return_hessian==False)&(return_gradient==True):
        return phi.astype(float), phi_i.astype(float)

    elif (return_hessian==True)&(return_gradient==True):
        return phi.astype(float), phi_i.astype(float), phi_ij.astype(float)







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


