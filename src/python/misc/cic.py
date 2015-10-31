import numpy as np





def cic(position, nbin):
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

    ''' ->> CIC <<- '''
    # ->> particle indices <<- #
    xc=np.floor(pos)
    idx=xc.astype(int)

    idx[np.nonzero(idx<0)]=idx[np.nonzero(idx<0)]+nbin
    idx[np.nonzero(idx>=nbin)]=idx[np.nonzero(idx>=nbin)]-nbin

    # ->> weighting function <<- #
    dx=np.fabs(pos-xc)
    tx=np.fabs(1.-dx)

    #print idx.shape, xc.shape, dx.shape, tx.shape


    ''' ->> density estimation <<- '''
    mass=1.
    d=np.zeros((nbin, nbin, nbin))

    #->> index <<- #
    idx_1=idx+1
    idx_1[np.nonzero(idx_1<0)]=idx_1[np.nonzero(idx_1<0)]+nbin
    idx_1[np.nonzero(idx_1>=nbin)]=idx_1[np.nonzero(idx_1>=nbin)]-nbin

    idx_x, idx_y, idx_z=idx[...,0], idx[...,1], idx[...,2]
    idx_x1,idx_y1, idx_z1=idx_1[...,0], idx_1[...,1], idx_1[...,2]

    #->> resume density <<- #
    d[idx_x,idx_y,idx_z]+=tx[idx_x,idx_y,idx_z,0]*tx[idx_x,idx_y,idx_z,1]*tx[idx_x,idx_y,idx_z,2]
    d[idx_x1,idx_y,idx_z]+=dx[idx_x,idx_y,idx_z,0]*tx[idx_x,idx_y,idx_z,1]*tx[idx_x,idx_y,idx_z,2]
    d[idx_x,idx_y1,idx_z]+=tx[idx_x,idx_y,idx_z,0]*dx[idx_x,idx_y,idx_z,1]*tx[idx_x,idx_y,idx_z,2]
    d[idx_x1,idx_y1,idx_z]+=dx[idx_x,idx_y,idx_z,0]*dx[idx_x,idx_y,idx_z,1]*tx[idx_x,idx_y,idx_z,2]

    d[idx_x,idx_y,idx_z1]+=tx[idx_x,idx_y,idx_z,0]*tx[idx_x,idx_y,idx_z,1]*dx[idx_x,idx_y,idx_z,2]
    d[idx_x1,idx_y,idx_z1]+=dx[idx_x,idx_y,idx_z,0]*tx[idx_x,idx_y,idx_z,1]*dx[idx_x,idx_y,idx_z,2]
    d[idx_x,idx_y1,idx_z1]+=tx[idx_x,idx_y,idx_z,0]*dx[idx_x,idx_y,idx_z,1]*dx[idx_x,idx_y,idx_z,2]
    d[idx_x1,idx_y1,idx_z1]+=dx[idx_x,idx_y,idx_z,0]*dx[idx_x,idx_y,idx_z,1]*dx[idx_x,idx_y,idx_z,2]

    return d*mass
