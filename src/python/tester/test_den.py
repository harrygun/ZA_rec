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
    'import_format':   'cita_simulation_highres',
    #'import_format':   'cita_simulation',
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
    #import_type='all'
    import_type='field'


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

    elif p.import_format=='cita_simulation_highres':
        p.nbin=576
	p.boxsize=1000.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'
	fn_write=droot_part+'0.000xv0.dat.displaced.npz'

	fn_field_write=droot_field+'0.000xv0.dat.den.npz'

        #dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)
        pos, v=mio.read_cita_simulation(fn_part, p.nbin)
	del v

	if True:
	    pl.plot(pos[:,:,100,1], pos[:,:,100,2], 'k.')
	    pl.show()


        for i in range(3):
            xmax[i], xmin[i] = np.max(pos_[...,i]), np.min(pos_[...,i])
            pos[...,i]=(pos[...,i]-xmin[i])*p.boxsize/(xmax[i]-xmin[i]) 


	# ->> estimate mass resolution <<- #
	print 'pos.shape', pos.shape, type(pos[0,0])
	print 'pos type:', type(pos[0,0])

        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass

    else:
        raise Exception






    #->> density contrast <<- #
    npt=pos.shape[0]
    print 'npart=', npt, 'nbin=', p.nbin, 
    d=mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos, redshift=p.redshift, pmass=p.particle_mass)

    np.savez(fn_field_write, d=d)





    
    # ->> The End <<- #
    p.finalize()
