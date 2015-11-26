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

import lagrangian.lpt_rec as lrec
import misc.io as mio
import fourier.potential as ptt
import misc.file_import as fimp
import misc.ps as ps
import misc.cic as mcic
import fourier.psxi as psxi











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
    'do_LPT_rec':        True,
    'import_reced_data':   False,
    #-------------------------------#
    'do_testing': False,
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

    # ->> imporot data <<- #
    if p.import_format=='MIP':
        p.nbin=256
	p.boxsize=32.
        droot_part='/mnt/scratch-lustre/xwang/data/velinv/MIP/particle/'
        droot_field='/mnt/scratch-lustre/xwang/data/velinv/MIP/256/raw/'
        fn='REAL_00-00-00'
        dd=fimp.import_MIP_data(p, droot_part, droot_field, fn, import_data_type=import_type)

    # ->>
    if p.import_format=='gadget_DTFE':
        p.nbin=256
	p.boxsize=100.
        droot_part='/mnt/scratch-lustre/xwang/data/velinv/cmpc/gadget/'
        droot_field='/mnt/scratch-lustre/xwang/data/velinv/cmpc/field/'

	fn_part=droot_part+'snap100Mpc256_z0'
	fn_field=droot_field+'snap100Mpc256_z0.den'

        dd=fimp.import_gadget_DTFE(p, fn_part, fn_field, import_data_type=import_type)

    if p.import_format=='cita_simulation':
        p.nbin=256
	p.boxsize=512.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'
	fn_write=droot_part+'0.000xv0.dat.displaced.npz'

        dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)

	# ->> estimate mass resolution <<- #

        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass




    ''' ->> Performing Lagrangian reconstruction <<- '''
    do_poisson=True
    if do_poisson==True:
        pos, delta =dd
	print 'data shape:', pos.shape, delta.shape, delta.min(), delta.max(), delta.mean()


        phi=ptt.Poisson3d(delta, boxsize=p.boxsize, smooth_R=0., smooth_type=p.smooth_type)
	print 'phi shape', phi.shape,  phi.min(), phi.max(), phi.mean()


        if p.save_data:

	    fn_write=droot_part+'0.000xv0.phi.npz'
	    try: fn_write
	    except: pass

	    np.savez(fn_write, phi=phi)


    elif p.import_reced_data==True:
        #->> 
	f=np.load(fn_write)
        d_rec, d_disp, d_shift = f['d_rec'], f['d_disp'], f['d_shift']




    
    






    
    # ->> The End <<- #
    p.finalize()
