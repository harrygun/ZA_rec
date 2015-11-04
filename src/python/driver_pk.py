import matplotlib
import matplotlib.pyplot as pl
import matplotlib.colors as colors

import numpy as np
import array as arr
import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.myplot as mpl

import misc.io as mio
import misc.cyth.cgalio as cg
import misc.cic as mcic
import misc.ps as ps


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
    import_type='all'


    if p.import_format=='cita_simulation':
        p.nbin=256
	p.boxsize=1000.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'
	fn_write=droot_part+'0.000xv0.dat.displaced.npz'
	fn_rec_write=droot_part+'0.000xv0.dat.displaced.den.npz'

        #dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)
	#pos= np.load(fn_write)['pos_displaced'].reshape(p.nbin**3, 3)

        di=np.load(fn_field)['d']
	dr=np.load(fn_rec_write)['d']

	print 'density shape:', di.shape, dr.shape


    # ->> power spectrum measurement <<- #
    ki, pki=ps.pk(di, boxsize=p.boxsize)
    kr, pkr=ps.pk(dr, boxsize=p.boxsize)

    print 'ki-kr:', ki-kr


    if True:
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
    
        ax[0].loglog(ki, pki, 'r--')
        ax[0].loglog(kr, pkr, 'k-')

        ax[1].plot(ki, pkr/pki) 
	ax[1].set_xscale("log")

        fig.savefig(root+'figure/ps_comp.png')
	pl.show()
