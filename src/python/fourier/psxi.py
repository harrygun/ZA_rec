import numpy as np
#import scipy.fftpack as sfft

import genscript.myarray as ar
import genscript.fft as fft





def pk(d, boxsize=1000.):
    # ->> obtain the angular averaged power spectrum <<- #

    # ->> FFT first <<- #
    dk = np.rfftn(d)
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

    print k2.shape

      


    return







def xi(d, boxsize=1000.):

    return
