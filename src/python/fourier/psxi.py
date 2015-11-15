import numpy as np
import scipy.fftpack as sfft

import genscript.myarray as ar
import genscript.fft as fft





def pk(d, boxsize=1000.):
    # ->> obtain the angular averaged power spectrum <<- #

    # ->> FFT first <<- #
    dk = np.fft.rfftn(d)
    s = d.shape
    sk = dk.shape

    dk2 = (dk*np.conjugate(dk)).astype(np.float32)

    # ->>  Fourier mode <<- #
    ng = d.shape[0]
    ndim=len(d.shape)
    kmin = 2.*np.pi/float(boxsize)

    # ->> Fourier space arguments <<- #
    kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng/2+1].astype(float)
    ki = 2.*np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), np.sin(kz*kmin/2.)])
    k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)
    k=np.sqrt(k2)

    print k2.shape, k.shape


    # ->> k bins <<- #
    bin2f=1/16.
    bedges=kmin*2.*np.arange(-bin2f/2., np.log(k[-1]/kmin)/np.log(2.), bin2f )


    #->> average over given k amplitude <<- #
    for i in range(len(bedges)-1):
        #->> 


    return







def xi(d, boxsize=1000.):

    return
