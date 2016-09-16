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
        if(R!=0):
            kr=k*R
            return 3.*(np.sin(kr)-kr*np.cos(kr) )/kr**3 
	else:
	    return 1.


def vsigma2_integrate(p, R=0., window_type='Gauss'):

    sf = lambda lgk: p.pk(10.**lgk)*wd(10.**lgk, R, window_type)**2*10.**(lgk)\
                        *np.log(10.)/2./np.pi**2


    s2=intg.quad(sf, -4., 3., limit=np.int(1e8) )

    return s2


def vsig2(p, z, R=0., window_type='Gauss'):
    s2=vsigma2_integrate(p, R=R, window_type=window_type)
    #return  (p.pk.D1(z)*p.pk.f(z)*p.dist.mH(z))**2*s2[0]
    return  s2[0]





param_dict={
    'power_spectrum_fname': '/Users/wangxin/workspace/code/pyport/default_data/pk/fiducial/fiducial_matterpower.dat',
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
    smooth_type= 'Tophat'  #'Tophat'

    vsig=np.sqrt(vsig2(p, z, R=R, window_type=smooth_type) )
    print 'sigma(z={0}, R={1})={2}'.format(z, R, vsig)

    root='../../workspace/'




    ''' ->> End of initialization <<- '''
    """
    nR=100
    Rlist=np.linspace(0, 20, nR)
    sigR=np.zeros(nR)


    for i in range(nR):
        R=Rlist[i]
        sigR[i]=np.sqrt(vsig2(p, 0., R=R, window_type='Tophat') )

    np.savetxt('sigmav_R.txt', (Rlist, sigR))
    """



    na=500
    alist=np.linspace(1e-2, 1, na)
    zlist=np.linspace(1e-2, 1, na)
    vsig_list=np.zeros(na)

    for i in range(na):
        z=1./alist[i]-1.
	zlist[i]=z
        vsig_list[i]=p.pk.D1(z)*p.pk.f(z)*p.dist.mH(z)*vsig

    np.savetxt('sigmav_Z.txt', zip(*(alist, vsig_list)))


    pl.plot(alist, vsig_list)
    pl.show()



    
    p.finalize()
