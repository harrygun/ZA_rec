import numpy as np
import pylab as pl
import sympy as sp
import gc, sys, os
import scipy.integrate as intg

import genscript.progcontrol as pc
from genscript.extendclass import *



def wd(k, R, window_type):
    if window_type=='Gauss':
        return np.exp(-(k*R)**2./2.)

    elif window_type=='Tophat':
        kr=k*R
        return 3.*(np.sin(kr)-kr*np.cos(kr) )/kr**3 


def vsigma2_integrate(p, R=0., window_type='Gauss'):

    #sf = lambda lgk: p.pk(10.**lgk)*wd(10.**lgk, R, window_type)**2*10.**(3.*lgk)\
    #                    *np.log(10.)/2./np.pi**2
    sf = lambda lgk: p.pk(10.**lgk)*wd(10.**lgk, R, window_type)**2*10.**(lgk)\
                        *np.log(10.)/2./np.pi**2

    s2=intg.quad(sf, -4., 3., limit=np.int(1e8) )

    return s2


def vsig2(p, z, R=0., window_type='Gauss'):
    s2=vsigma2_integrate(p, R=R, window_type=window_type)
    return  (p.pk.D1(z)*p.pk.f(z)*p.dist.mH(z))**2*s2[0]





param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'cosmology_parameter_fname': 'parameters/cosparameter.cfg',
    'cosmology_parameter_sec': 'Cosmology_Parameters',
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
    ''' ->> Initialization <<- '''
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    # ->> sigma estimation <<- #
    z=0.  # redshift
    R=0.  # Mpc/h
    smooth_type= 'Gauss'  #'Tophat'

    vsig=np.sqrt(vsig2(p, z, R=R, window_type=smooth_type) )
    print 'sigma(z={0}, R={1})={2}'.format(z, R, vsig)

    root='../../workspace/'

    ''' ->> End of initialization <<- '''

    
    p.finalize()
