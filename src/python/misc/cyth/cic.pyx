cimport numpy as np
import numpy as np



cpdef cic_sum(int npart, int nbin, np.ndarray[np.double_t, ndim=3]d, \
              np.ndarray[np.double_t, ndim=2]pos):

    cdef:
        int i, j, k, i1, j1, k1, ip
        double dx, dy, dz, tx, ty, tz


    # ->> CIC sum over all particles
    for ip in range(npart):
        i, j, k=<int>pos[ip,0], <int>pos[ip,1], <int>pos[ip,2]
        i1, j1, k1=i+1, j+1, k+1

        # ->> perodic condition <<- #
        if i<0: i=i+nbin
        if j<0: j=j+nbin
        if k<0: k=k+nbin
      
        if i>=nbin: i=i-nbin
        if j>=nbin: j=j-nbin
        if k>=nbin: k=k-nbin

        if i1<0: i1=i1+nbin
        if j1<0: j1=j1+nbin
        if k1<0: k1=k1+nbin
      
        if i1>=nbin: i1=i1-nbin
        if j1>=nbin: j1=j1-nbin
        if k1>=nbin: k1=k1-nbin

        #->> kernel <<- #
        dx=np.fabs(pos[ip,0]-<double>i)
        dy=np.fabs(pos[ip,1]-<double>j)
        dz=np.fabs(pos[ip,2]-<double>k)
	
        tx=np.fabs(1.0-dx)
        ty=np.fabs(1.0-dy)
        tz=np.fabs(1.0-dz)

        #->> summation <<- #
        d[i,j,k]+=tx*ty*tz
        d[i1,j,k]+dx*ty*tz
        d[i,j1,k]+tx*dy*tz
        d[i1,j1,k]+=dx*dy*tz
        d[i,j,k1]+=tx*ty*dz
        d[i1,j,k1]+=dx*ty*dz
        d[i,j1,k1]+=tx*dy*dz
        d[i1,j1,k1]+=dx*dy*dz

    return






cpdef density(long long npart, int nbin, np.ndarray[np.double_t, ndim=2]pos, double mass): 

    cdef:
        int i
        long long Ngrid[3]
        double masstot

    # ->> initialization <<- #
    delta=np.zeros((nbin, nbin, nbin))

    for i in range(3):
        Ngrid[i]=<long long>nbin 

    # ->> call density <<- #
    masstot=density(<float *>pos.data, <float *>delta.data, <long long> npart, \
                    <long long> Ngrid[3], <double> mass)

    # ->>
    '''
    aa=0;
    x1=(xmax-xmin)/Ngridx/HubbleParam
    y1=(ymax-ymin)/Ngridy/HubbleParam
    z1=(zmax-zmin)/Ngridz/HubbleParam
    for(i=0;i<=Ngridx-1;i++)
       for(j=0;j<=Ngridy-1;j++)
          for(k=0;k<=Ngridz-1;k++){
             delta[i][j][k]=delta[i][j][k]/(x1*y1*z1)/rhom*1e10/HubbleParam-1.
             aa+=delta[i][j][k]
             }
    printf("mean density %lg\n",aa/Ngridx/Ngridy/Ngridz);
    '''

    return delta
