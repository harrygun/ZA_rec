import matplotlib
#matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as pl
import matplotlib.colors as colors

import os
import numpy as np
import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.read as read
import genscript.myplot as mpl
import genscript.read as rd

import misc.io as mio
import fourier.potential as ptt
import misc.file_import as fimp
import misc.ps as ps
import misc.cic as mcic

import fourier.psxi as psxi 
import misc.power as pk_cor










param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 15.,
    'smooth_type': 'Gaussian', 
    'smooth_R_list_type':  'linear', 
    'boxsize': 32.,
    'nbin':    256, 
    'particle_mass':    1.e5,
    #'import_format':   'gadget_DTFE',
    'import_format':   'cita_simulation',
    'save_data':    True,
    'redshift':     0.,
    }

prog_control={
    #-------------------------------#
    'do_pk':        True,
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/result/'


    ''' -------------------------------------------------
	-------------------------------------------------     
    '''
    import_type='all'
    #import_type='field'


    if p.import_format=='cita_simulation':
        p.nbin=256
	p.boxsize=1000.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'
	fn_write=droot_part+'0.000xv0.dat.displaced.npz'

        dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)

	# ->> estimate mass resolution <<- #

        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass

    else:
        raise Exception



    #->> density contrast <<- #
    delta =dd[1]


    ''' ->> analysis the data <<- '''
    do_powerspectrum = False
    do_correlation_function = True

    if do_powerspectrum==True:
        print 'boxsize=', p.boxsize
    
	k, pk=psxi.pk(delta, boxsize=p.boxsize)
        k_ori, pk_ori=ps.pk(delta, boxsize=p.boxsize)

	print 'k/k_ori shape:', k.shape, k_ori.shape


        if True:
            nplt, ncol = 1, 1
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)

            ax[0].loglog(k_ori, pk_ori, 'r--')
            ax[0].loglog(k, pk, 'k-')

            #ax[1].plot(k, pk_ori/pk, 'k-') 
            #ax[1].set_xscale("log")
    
            pl.show()


    if do_correlation_function==True:

        #r, xi=psxi.xi(delta, boxsize=p.boxsize)
        #xi=psxi.xi(delta, boxsize=p.boxsize)

	r, xi =pk_cor.corfunk(delta, boxsize=p.boxsize )


	pl.plot(r, r**2.*xi)
	pl.show()






    
    # ->> The End <<- #
    p.finalize()
