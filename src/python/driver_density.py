import matplotlib
#matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
import matplotlib.colors as colors

import numpy as np
import array as arr
import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *

import misc.io as mio
import misc.cyth.cgalio as cg
import misc.cic as mcic


import misc.file_import as fimp






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 15.,
    'smooth_type': 'Gaussian', 
    'smooth_R_list_type':  'linear', 
    'boxsize': 32.,
    'nbin':    256, 
    #'import_format':   'gadget_DTFE',
    'import_format':   'cita_simulation',
    'save_displaced_particles':    True,
    'particle_mass':    1.e5,
    'redshift':     0.,
    }

prog_control={
    #-------------------------------#
    'do_LPT_rec': True,
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
    #import_type='field'
    import_type='particle'

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
	p.boxsize=1000.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'
	fn_write=droot_part+'0.000xv0.dat.displaced.npz'
	fn_rec_write=droot_part+'0.000xv0.dat.displaced.den.npz'

	fn_field_write=droot_field+'0.000xv0.dat.den.npz'

        dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)
	#pos= np.load(fn_write)['pos_displaced'].reshape(p.nbin**3, 3)

        pos=np.copy(np.swapaxes(dd[0].reshape(3, p.nbin**3), 0, 1))
	print 'pos.shape', pos.shape, type(pos[0,0])

        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass


    '''----------------------------------------------------------------------------- '''


    npt=pos.shape[0]
    print 'npart=', npt, 'nbin=', p.nbin, 

    d=mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos, redshift=p.redshift, pmass=p.particle_mass)
    print 'd shape:', d.shape

    print 'd min/max:', d.min(), d.max(), d.mean()

    np.savez(fn_field_write, d=d)


    if True:
        fig=pl.figure(figsize=(20, 20))
        ax=fig.add_subplot(111)
    
        data=d[...,100]+1.001
        ax.imshow(np.flipud(data), norm=colors.LogNorm(vmin=data.min(),vmax=data.max()) )
    
        #fig.savefig('cita_test.png')
	pl.show()
