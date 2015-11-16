import numpy as np
import scipy.fftpack as sfft

import genscript.myarray as ar
import genscript.fft as fft
import genscript.myarray as mar





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
    #kx, ky, kz=np.mgrid[0:ng, 0:ng, 0:ng/2+1].astype(float)
    #ki = 2.*np.array([np.sin(kx*kmin/2.), np.sin(ky*kmin/2.), np.sin(kz*kmin/2.)])
    #k2=4.*(np.sin(kx*kmin/2.)**2.+np.sin(ky*kmin/2.)**2.+np.sin(kz*kmin/2.)**2.)

    #ki = np.array([kx*kmin, ky*kmin, kz*kmin])
    #k2=(kx*kmin)**2.+(ky*kmin)**2.+(kz*kmin)**2.


    ki_list = [ sfft.fftfreq(s[i],d=1./float(s[i]))*kmin for i in range(2) ]\
              + [sfft.rfftfreq(s[-1], d=1./float(s[-1]))*kmin ]
    print 'kx, ky, kz max:', np.max(ki_list[0]), np.max(ki_list[1]), np.max(ki_list[2])

    ki=mar.meshgrid(*ki_list)
    k2=ki[0]**2.+ki[1]**2.+ki[2]**2.



    _k=np.sqrt(k2).flatten()
    kmax=_k.max()
    index = np.argsort(_k)

    k=_k[index]
    dk2=_dk2[index]

    print 'kmin/kmax=', kmin, kmax, k[-1]
    print 'k2.shape=', k2.shape, 'k.shape=', k.shape



    # ->> k bins <<- #
    bin2f=1./16.
    bedges=kmin*2.**np.arange(-bin2f/2., np.log(kmax/kmin)/np.log(2.), bin2f)
    print 'k-space edges', len(bedges), kmin, bin2f, kmax, k[-1]


    #->> cuts <<- #
    cuts = np.searchsorted(k,bedges)
    #print 'error: ', (k[cuts]-bedges)/bedges


    # ->> 
    numinbin = np.zeros(bedges.shape)
    pk = np.zeros(bedges.shape)
    kmean = np.zeros(bedges.shape)
    nbins = len(bedges)

    _kz0=np.ones(ki[-1].flatten().shape)
    _kz0[np.where(ki[-1].flatten()==0.)]-=0.5
    kz0=_kz0[index]


    #->> average over given k amplitude <<- #
    for i in range(len(bedges)-1):
        #->> 
        if (cuts[i+1]>cuts[i]):
            numinbin[i] = np.sum(kz0[cuts[i]:cuts[i+1]])
            pk[i] = np.sum(kz0[cuts[i]:cuts[i+1]]*dk2[cuts[i]:cuts[i+1]])
            kmean[i] = np.sum(kz0[cuts[i]:cuts[i+1]]*k[cuts[i]:cuts[i+1]])

    #print 'numinbin:', numinbin, len(numinbin)
    #quit()

    wn0 = np.where(numinbin > 0.)[0]
    pk = pk[wn0]; kmean = kmean[wn0]; numinbin=numinbin[wn0]
    pk /= numinbin
    kmean /= numinbin

    #print 'wn0:', wn0

    #quit()

    pk *= boxsize**3/np.prod(np.array(s).astype(float))**2

    return kmean, pk







def xi(d, k=None, ps=None, boxsize=1000.):

    #->> do power spectrum first <<- #
    if ((k==None)|(ps==None)):
        k, ps=pk(d, boxsize=boxsize)

    #->> 

    

    return


