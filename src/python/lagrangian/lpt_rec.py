import os
import numpy as np
import pylab as pl
#import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.interpolation as itp

import fourier.potential as ptt




def get_displacement(p, m, boxsize, smooth_R=None, smooth_type=None):
    # ->> 
    phi=ptt.Poisson3d(m, boxsize=boxsize, smooth_R=smooth_R, smooth_type=smooth_type)
    si = np.array(np.gradient(phi))

    # ->> return the inverse of ZA <<- #
    return phi, si




def lag_rec(p, mpart, map, boxsize, smooth_R=None, smooth_type=None): 
    # -> reconstruction <<- #

    phi, si=get_displacement(p, map, boxsize, smooth_R=smooth_R, smooth_type=smooth_type)

    # ->> setup interpolation of si <<- #
    x, y, z=
    data =[np.array([x, y, z, si[i]]) for i in range(3)]

    sitp=[]
    for i in range(3):
        sitp.append(itp.griddata_interpolater(data[i], 3, method='cubic', fill_value=0))

    # ->> move particles <<- #
    print '->> shape:', phi.shape, si.shape


    return phi, si
