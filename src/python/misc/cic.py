import numpy as np
import cyth.cic as cicyth



_cic_type_ = 'C version'
#_cic_type_ = 'Cython version'


def cic(cp, npart, nbin, boxsize, position, pmass=1e9):
    ''' ->> return the cloud-in-cell density <<- '''

    # ->> regulate particle positions <<- #
    xmin, xmax =np.zeros(3), np.zeros(3)
    dl=np.zeros(3)
    pos=np.zeros(position.shape)

    for i in range(3):
        xmax[i], xmin[i] = np.max(position[...,i]), np.min(position[...,i])
        #xmax, xmin = np.max(xmax_, boxsize), np.min(xmin, )

        dl[i]=(xmax[i]-xmin[i])/float(nbin)
        pos[...,i]= (position[...,i]-xmin[i])/dl[i]

    #print 'dl=', dl

    ''' ->> CIC density estimation <<- '''

    # ->> summation over all particles <<- #
    d=np.zeros((nbin, nbin, nbin))

    if _cic_type_ == 'Cython version':

        d=np.zeros((nbin, nbin, nbin))
        cicyth.cic_sum(npart, nbin, d, pos)

        #print npart, nbin, d.shape
        #print 'nonzero:', len(np.where(d!=0)[0])
        #print d[np.where(d!=0)]


    elif _cic_type_ == 'C version':

        d=cicyth.density_cyth(npart, nbin, pos, pmass)

        rhom=(2.7752e11)*cp.h**2.*cp.omem 
        x1=np.array([(xmax[i]-xmin[i])/np.float(nbin)/cp.h  for i in range(3) ])

        d=d/np.prod(x1)/cp.h**2./rhom*1e10/cp.h-1.

        print 'get density:', d.shape, rhom



    return d
