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
                                 np.ndarray[np.float32_t, ndim=3, mode='c']shifted, \
		                 np.ndarray[np.float32_t, ndim=2, mode='c']grid  ):
    ''' ->> pos[npart, 3], grid[3, ngrid], si[3, ngrid, ngrid, ngrid] <<- '''

    cdef:
        int n, ip, i, j, k, i1, j1, k1, ni, nj, nk
        float dx, dy, dz, ddx, ddy, ddz, dp[3], x[3], v[2][2][2], dsi[3]

    # ->> uniform grid <<- #
    dp[0],dp[1],dp[2]=grid[0,1]-grid[0,0], grid[1,1]-grid[1,0], grid[2,1]-grid[2,0]


    for ip in range(npart):
        #->>  
        i=np.searchsorted(grid[0], pos[ip,0])
        j=np.searchsorted(grid[1], pos[ip,1])
        k=np.searchsorted(grid[2], pos[ip,2])

        dx=(pos[ip,0]-grid[0,i])/ddx
        dy=(pos[ip,1]-grid[1,j])/ddy
        dz=(pos[ip,2]-grid[2,k])/ddz


        #->>
	for n in range(3):

            for ni in range(2):
                i1=i+1
                if i1>ngrid

                for nj in range(2):
                    for nk in range(2):
                        v[ni][nj][nk]=si[n,i+1,j+1,k+1]

	    dsi[n]=trilinear(pos[3], float dp[3], float v[3][3][3])



        shifted[ip,n]=pos[ip,n]+si[]
         

    return




