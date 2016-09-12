import numpy as np
import pylab as pl
import sympy as sp
import gc, sys, os




def wd(k, R, window_type):
    if window_type=='Gauss':
        return np.exp(-(k*R)**2./2.)

    elif window_type=='Tophat':
        kr=k*R
        return 3.*(np.sin(kr)-kr*np.cos(kr) )/kr**3 


def vsigma2_integrate(p, R=0., window_type='Gauss'):

    sf = lambda lgk: p.pk(10.**lgk)*wd(10.**lgk, R, window_type)**2*10.**(3.*lgk)\
                        *np.log(10.)/2./np.pi**2

    s2=intg.quad(sf, -3., 3., limit=np.int(1e8) )

    return s2


def vsig2(p, z, R=0., window_type='Gauss'):
    s2=vsigma2_integrate(p, R=R, window_type=window_type)
    return  p.pk.D1(z)**2*s2[0]






if __name__=='__main__':
    ''' ->> Initialization <<- '''
    sec=''
    init_dict=myDict({})
    p=pc.prog_init(section=sec, **init_dict)

    # ->> sigma estimation <<- #
    z=0.  # redshift
    R=0.  # Mpc/h
    smooth_type= 'Gauss'  #'Tophat'

    vsig=np.sqrt(vsig2(p, z, R=R, window_type=p.smooth_type) )
    print 'sigma(z={0}, R={1})={2}'.format(z, R, vsig)

    root='../../workspace/'

    ''' ->> End of initialization <<- '''

    
    p.finalize()
