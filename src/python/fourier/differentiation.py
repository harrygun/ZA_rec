import numpy as np
import scipy.fftpack as sfft
import scipy.ndimage.filters as sft

from genscript.extendclass import *
import genscript.myarray as ar




def diff1(d):
    # ->> 1st order differentiation using FFT


    return


def vfield_gradient_real(v, nbin, dl=1.):
    dvij=[ np.gradient(v[:,:,:,i] ) for i in range(3) ]    

    dvij=np.array(dvij)/dl
    print 'gradient shape: {0}'.format(dvij.shape)

    return dvij


def vfield_gradient_fft(v, nbin):
    s=v.shape
    print 'v shape:', s

    ng = s[1]
    kmin = 2.*np.pi/float(ng)

    # ->> Do NOT shift freq for `sin(k)' kernel
    kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng/2+1].astype(float)
    ik=2j*np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), \
                    np.sin(kz*kmin/2.)])

    
    #dk=np.array([np.fft.rfftn(v[...,i]) for i in range(3)])
    dk=np.array([np.fft.rfftn(v[...,i]) for i in range(3)])

    gk=np.array([ [ik[i]*dk[j] for i in range(3) ] for j in range(3)])
    print 'gk shape:', gk.shape

    gk[:,:,0,0,0]=0

    dvij=np.array([ [np.fft.irfftn(gk[i,j]) for i in range(3) ] for j in range(3)] )
    #dvij=np.fft.irfftn(gk, axes=[2,3,4]) 
    print dvij.shape

    return dvij


def vorticity(v, boxsize, nbin):
    dl=boxsize/nbin

    dv=vfield_gradient_real(v, nbin, dl=dl)
    #dv=vfield_gradient_fft(v, nbin)

    wx=dv[1][2]-dv[2][1]
    wy=dv[2][0]-dv[0][2]
    wz=dv[0][1]-dv[1][0]

    ww=np.sqrt(wx**2.+wy**2.+wz**2.)

    return np.array([wx, wy, wz]), ww
