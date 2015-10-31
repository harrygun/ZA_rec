import numpy as np
import cyth.cic as ccic




def cic(npart, nbin, position, pmass=1e9):
    ''' ->> return the cloud-in-cell density <<- '''

    # ->> regulate particle positions <<- #
    xmin, xmax =np.zeros(3), np.zeros(3)
    dl=np.zeros(3)
    pos=np.zeros(position.shape)

    for i in range(3):
        xmax[i], xmin[i] = np.max(position[...,i]), np.min(position[...,i])
        dl[i]=(xmax[i]-xmin[i])/float(nbin)
        pos[...,i]= (position[...,i]-xmin[i])/dl[i]

    #print 'dl=', dl

    ''' ->> CIC density estimation <<- '''
    mass=pmass
    d=np.zeros((nbin, nbin, nbin))

    # ->> summation over all particles <<- #
    ccic.cic_sum(npart, nbin, d, pos)
    #print npart, nbin, d.shape

    #print 'nonzero:', len(np.where(d!=0)[0])
    #print d[np.where(d!=0)]


    return d*mass
