import os
import numpy as np
import pylab as pl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power




def Testing_portal(p):

    return




param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 0,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    'do_testing': False, 
    #-------------------------------#
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/result/'



    # ->> smoothing list <<- #
    p.z= 0.

    # ->> Rmin, Rmax:  only integer numbers for convinience of file name <<- #



    ''' -------------------------------------------------
               ->>      Making some plots <<- 
	-------------------------------------------------     
    '''






    ''' -------------------------------------------------
               ->>      do some testing      <<- 
	-------------------------------------------------     
    '''
    if p.do_testing==True:
        Testing_portal(p)



    
    # ->> The End <<- #
    p.finalize()
