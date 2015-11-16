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

    _dk2 = (dk*np.conjugate(dk)).flatten()  #.astype(np.float32)

    # ->>  Fourier mode <<- #
    ng = d.shape[0]
    ndim=len(d.shape)
    kmin = 2.*np.pi/float(boxsize)

    # ->> Fourier space arguments <<- #
    kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng/2+1].astype(float)
    ki = 2.*np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), np.sin(kz*kmin/2.)])
    k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)

    _k=np.sqrt(k2).flatten()
    kmax=_k.max()
    index = np.argsort(_k)

    k=_k[index]
    dk2=_dk2[index]

    print 'kmin/kmax=', kmin, kmax
    print 'k2.shape=', k2.shape, 'k.shape=', k.shape


    # ->> k bins <<- #
    bin2f=1/16.
    bedges=kmin*2.**np.arange(-bin2f/2., np.log(kmax/kmin)/np.log(2.), bin2f)
    print 'k-space edges', bedges


    #->> cuts <<- #
    cuts = np.searchsorted(k,bedges)
    #print 'error: ', (k[cuts]-bedges)/bedges


    # ->> 
    numinbin = np.zeros(bedges.shape)
    pk = np.zeros(bedges.shape)
    kmean = np.zeros(bedges.shape)
    nbins = len(bedges)


    _kz0=np.ones(ki[-1].shape)
    _kz0[np.where(ki[-1].flatten()==0.)]-=0.5
    kz0=_kz0[index]


    #->> average over given k amplitude <<- #
    for i in range(len(bedges)-1):
        #->> 
        if (cuts[i+1]>cuts[i]):
            numinbin[i] = np.sum(kz0[cuts[i]:cuts[i+1]])
            pk[i] = np.sum(kz0[cuts[i]:cuts[i+1]]*dk2[cuts[i]:cuts[i+1]])
            kmean[i] = np.sum(kz0[cuts[i]:cuts[i+1]]*k[cuts[i]:cuts[i+1]])
	else:
	    print 'cut[i+1]<=cut[i]', i, cuts[i+1], cuts[i]


    wn0 = np.where(numinbin > 0.)[0]
    pk = pk[wn0]; kmean = kmean[wn0]; numinbin=numinbin[wn0]
    pk /= numinbin
    kmean /= numinbin

    return kmean, pk







def xi(d, boxsize=1000.):

    return


