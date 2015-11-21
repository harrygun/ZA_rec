cimport numpy as np
import numpy as np






cdef float trilinear(float pos[3], float dp[3], float v[2][2][2]):

    cdef:
        float c00, c01, c10, c11, c0, c1, c

    c00=v[0][0][0]*(1.-dp[0])+v[1][0][0]*dp[0]
    c10=v[0][1][0]*(1.-dp[0])+v[1][1][0]*dp[0]
    c01=v[0][0][1]*(1.-dp[0])+v[1][0][1]*dp[0]
    c11=v[0][1][1]*(1.-dp[0])+v[1][1][1]*dp[0]
    
    c0=c00*(1.-dp[1]) + c10*dp[1]
    c1=c01*(1.-dp[1]) + c11*dp[1]
    
    c=c0*(1.-dp[2])+c1*dp[2]

    return c




cdef pmove(int npart, int ngrid, np.ndarray[np.float32_t, ndim=2, mode='c']pos, \
                                 np.ndarray[np.float32_t, ndim=4, mode='c']si,  \
                                 np.ndarray[np.float32_t, ndim=2, mode='c']shifted, \
		                 np.ndarray[np.float32_t, ndim=2, mode='c']grid  ):
    ''' ->> pos[npart, 3], grid[3, ngrid], si[3, ngrid, ngrid, ngrid] <<- '''

    cdef:
        int n, ip, i, j, k, idx[3], i1, j1, k1, ni, nj, nk
        float delta_grid[3], dp[3], x[3], v[2][2][2], dsi

    print 'grid max:', grid[:,-1], grid.max()


    for ip in range(npart):
        #->>  
        i=np.searchsorted(grid[0], pos[ip,0])
        j=np.searchsorted(grid[1], pos[ip,1])
        k=np.searchsorted(grid[2], pos[ip,2])

        idx[0], idx[1], idx[2]=i, j, k

        print 'ip:', ip, i, j, k

        #->> preparing the interpolation <<- #
        for n in range(3):
            # ->> uniform grid <<- #
            delta_grid[n] = grid[n,1]-grid[n,0]
	    #->> particle position <<- #
            x[n]=pos[ip,n]

            print 'n:', n, 'idx[n]', idx[n]

            #->> distance to the grid <<- #
            dp[n]=(x[n]-grid[n,idx[n]])/delta_grid[n]


        print '\n'

        # ->> set vertices <<- #
        for n in range(3):
            for ni in range(2):
                i1=i+ni
                if i1>ngrid-1: i1=i1-ngrid

                for nj in range(2):
                    j1=j+ni
                    if j1>ngrid-1: j1=j1-ngrid

                    for nk in range(2):
                        k1=k+ni
                        if k1>ngrid-1: k1=k1-ngrid

                        v[ni][nj][nk]=si[n,k1,j1,k1]

            # ->> interpolation <<- #
            dsi=trilinear(x, dp, v)


            shifted[ip,n]=pos[ip,n]+dsi
         

    return




cpdef particle_move(p, pos, si):

    _grid_=np.linspace(0, p.boxsize, p.nbin, endpoint=False).astype(np.float32)
    grid=np.ascontiguousarray(np.vstack( (_grid_, _grid_, _grid_) ), dtype=np.float32)
    print 'grid shape:', grid.shape

    shifted=np.ascontiguousarray(np.zeros(pos.shape), dtype=np.float32)
    print 'cython partmoving shape:', shifted.shape, grid.shape, pos.shape, si.shape

    # ->> call cdef routine <<- #
    pmove(<int>p.nbin**3, <int>p.nbin, np.ascontiguousarray(pos, dtype=np.float32),\
          np.ascontiguousarray(si, dtype=np.float32), shifted, grid)

    return shifted
