import os
import numpy as np
import pylab as pl
#import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.interpolation as itp
import genscript.myutlis as mut
import genscript.myarray as mar

import fourier.potential as ptt




def get_displacement(p, m, boxsize, smooth_R=None, smooth_type=None):
    # ->> 
    phi=ptt.Poisson3d(m, boxsize=boxsize, smooth_R=smooth_R, smooth_type=smooth_type)

    # ->> inverse of ZA: Psi_i = -partial_i * Laplacian^{-1} delta(tau) <<- #
    dl=p.boxsize/np.float(p.nbin)
    si = np.array(np.gradient(phi))/dl

    print 'dl=', dl, '|si|<', np.max(np.fabs(si))

    return phi, si




def lag_rec(p, mpart, map, smooth_R=None, smooth_type=None): 
    # -> reconstruction <<- #

    phi, si=get_displacement(p, map, p.boxsize, smooth_R=smooth_R, smooth_type=smooth_type)
    print '->> phi, si shape:', phi.shape, si.shape

    # ->>
    moved=np.zeros((3, p.nbin, p.nbin, p.nbin))
    displace_interpolation=False

    if displace_interpolation==False:

        for i in range(3):
            dis_=mpart[i]+si[i]

            # ->> periodic boundary <<- #
	    dis_[np.where(dis_<0.)] = p.boxsize+dis_[np.where(dis_<0.)]
	    dis_[np.where(dis_>p.boxsize)] = dis_[np.where(dis_>p.boxsize)]-p.boxsize
	    moved[i] = np.copy(dis_)

    	    print 'axis-', i, 'is done.'


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
	    moved[i] = np.copy(dis_)

    	    print 'axis-', i, 'is done.'


    return moved
