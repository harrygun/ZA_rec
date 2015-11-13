import os
import numpy as np
#import pylab as pl
import matplotlib.pyplot as pl
#import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.interpolation as itp
import genscript.myutlis as mut
import genscript.myarray as mar

import fourier.potential as ptt
import misc.cic as mcic



def get_ZA_displacement(p, m, smooth_R=None, smooth_type=None):
    # ->> Poisson solver <<- #

    phi=ptt.Poisson3d(m, boxsize=p.boxsize, smooth_R=smooth_R, smooth_type=smooth_type)

    if False:
        # ->> making plots <<- #
        pl.imshow(phi[:,:,100])
        pl.show()
        quit()

    # ->> inverse of ZA: Psi_i = -partial_i * Laplacian^{-1} delta(tau) <<- #
    dl=p.boxsize/np.float(p.nbin)
    si = np.array(np.gradient(phi))/dl

    print 'dl=', dl, '|si|<', np.max(np.fabs(si))

    if False:
        pl.hist(si.flatten(), bins=20)
        pl.show()
        #print 'histogram:', hist


    return phi, si



def shifted_ZA(p, si):
    ''' ->> shift uniform grid after obtaining '''

    lgrid=np.linspace(0, p.boxsize, p.nbin, endpoint=False)
    print 'lgrid shape:', lgrid.shape

    # ->>  <<- #
    grid=mar.meshgrid(lgrid, lgrid, lgrid)[::-1,...]
    print 'meshgrid shape:', grid.shape, si.shape


    shifted=np.zeros((3, p.nbin, p.nbin, p.nbin))
    for i in range(3):
        shift_=grid[i]+si[i]
    
        # ->> periodic boundary <<- #
        shift_[np.where(shift_<0.)] = p.boxsize+shift_[np.where(shift_<0.)]
        shift_[np.where(shift_>p.boxsize)] = shift_[np.where(shift_>p.boxsize)]-p.boxsize
        shifted[i]=np.copy(shift_)
    
        print 'shifting axis-', i, 'is done.'

    return shifted




def displaced_ZA(p, si, mpart):
    # -> reconstruction <<- #

    displaced=np.zeros((3, p.nbin, p.nbin, p.nbin))

    displace_interpolation=False
    if displace_interpolation==False:

        for i in range(3):
            dis_=mpart[i]+si[i]

            # ->> periodic boundary <<- #
	    dis_[np.where(dis_<0.)] = p.boxsize+dis_[np.where(dis_<0.)]
	    dis_[np.where(dis_>p.boxsize)] = dis_[np.where(dis_>p.boxsize)]-p.boxsize
	    displaced[i] = np.copy(dis_)

    	    print 'displacing axis-', i, 'is done.'


    elif displace_interpolation==True:

        # ->> setup interpolation of si <<- #
        x_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
        y_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
        z_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
    
        x, y, z=mar.meshgrid(x_, y_, z_)
        data =[np.array([x, y, z, si[i]]) for i in range(3)]
    
        print 'interpolation data shape:', x_.shape, x.shape, data[0].shape
    
        sitp=[]
        for i in range(3):
            sitp.append(itp.griddata_interpolater(data[i], 3, method='linear', fill_value=0))
    
    
        # ->> move particles <<- #
        for i in range(3):
            dis_=mpart[i]+sitp[i]( (mpart[0], mpart[1], mpart[2]) )

            # ->> periodic boundary <<- #
	    dis_[np.where(dis_<0.)] = p.boxsize+dis_[np.where(dis_<0.)]
	    dis_[np.where(dis_>p.boxsize)] = dis_[np.where(dis_>p.boxsize)]-p.boxsize
	    displaced[i] = np.copy(dis_)

    	    print 'axis-', i, 'is done.'

    return displaced





''' ->> Density Reconstruction Wrapper <<- '''
_rect_type_list_=['ZA_displaced', 'ZA_displaced_shifted']

def lag_rec_ZA(p, mpart, dmap, smooth_R=None, smooth_type=None, rect_type='ZA_displaced'): 
    # -> reconstruction of initial state, return density map <<- #

    # ->> get displacement field <<- #
    phi, si=get_ZA_displacement(p, dmap, smooth_R=smooth_R, smooth_type=smooth_type)
    print '->> phi, si shape:', phi.shape, si.shape

    npt=p.nbin**3.

    if rect_type=='ZA_displaced_shifted':
        #->> get particles <<- #
        pos_disp = np.copy(np.swapaxes(displaced_ZA(p, si, mpart).reshape(3, p.nbin**3), 0, 1))
	pos_shift = np.copy(np.swapaxes(shifted_ZA(p, si).reshape(3, p.nbin**3), 0, 1))

        _test_draw_ = True
	if _test_draw_:
            #pl.plot( )
	    pass



	print 'particle shape:', mpart.shape, pos_disp.shape, pos_shift.shape

        #->> converting density map <<- #
	d_disp = mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos_disp, pmass=p.particle_mass)
	d_shift= mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos_shift, pmass=p.particle_mass)

	d_rec=d_disp-d_shift

        return d_rec, d_disp, d_shift, pos_disp, pos_shift


    elif rect_type=='ZA_displaced':

        return 


    else:
        raise Exception









