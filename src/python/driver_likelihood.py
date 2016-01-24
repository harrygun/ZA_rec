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
import genscript.myarray as mar

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
    'save_data':    False,
    'redshift':     0.,
    }

prog_control={
    #-------------------------------#
    'do_likelihood_rec':        True,
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
    # ->> imporot data <<- #
    import_type='particle'

    if p.import_format=='cita_simulation':
        p.nbin=256
	p.boxsize=512.
        droot_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
        droot_field='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'

	fn_part=droot_part+'0.000xv0.dat'
	fn_field=droot_field+'0.000xv0.dat.den.npz'

        dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type, renormalize_pos=False)

	# ->> estimate mass resolution <<- #
        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass

	# ->> import the initial density field <<- #
	fn_part_init=droot_part+'100.000xv0.dat'
        pos_init, v_init=mio.read_cita_simulation(fn_part_init, p.nbin)
	print 'init pos, v:', pos_init.shape, v_init.shape, pos_init.min(), pos_init.max()
         
    else:
        raise Exception



    ''' ->> Performing likelihood reconstruction <<- '''
    if p.do_likelihood_rec==True:

        # ->> obtain displacement field first <<- #
        pos=dd[0].reshape(p.nbin,p.nbin,p.nbin,3)
	pos_init=pos_init.reshape(pos.shape)
	disp=pos-pos_init

	print 'pos/pos_init/disp shape:', pos.shape, pos_init.shape, disp.shape
	print 'disp: ', disp.min(), disp.max()

        #->> another way of calculating displacement field <<- #
        lgrid=np.linspace(0, p.boxsize, p.nbin, endpoint=False)
        grid=np.rollaxis(mar.meshgrid(lgrid, lgrid, lgrid)[::-1,...], 0, 4)
        print 'lgrid shape:', lgrid.shape, 'meshgrid shape:', grid.shape

        err=pos_init-grid
	print 'err: ', err.max(), err.min()

	disp_uni=pos-grid

        if False:
	    nplots=2
            ax=mpl.mysubplots(nplots, ncol_max=2, subp_size=5)
	    n_bin=500

            ax[0].hist(disp[...,0].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='g')
            ax[0].hist(disp[...,1].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='r')
            ax[0].hist(disp[...,2].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='b')


            ax[1].hist(disp_uni[...,0].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='g')
            ax[1].hist(disp_uni[...,1].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='r')
            ax[1].hist(disp_uni[...,2].flatten(), bins=n_bin, range=[-20,20], histtype='step', color='b')


            pl.show()

	if True:
            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=1, 20
            
            ax[0].imshow(disp[axis,:,:,nsl]); 
            ax[1].imshow(disp_uni[axis,:,:,nsl]); 
            
            pl.tight_layout()
            pl.show()

        quit()

        # ->>  obtain ZA displacement <<- #
    """
        rect_type='ZA_displaced_shifted'

        _dd_ = lrec.lag_rec_ZA(p, pos, delta, smooth_R=p.smooth_R, smooth_type=p.smooth_type, rect_type=rect_type)

        # ->> unwrap density <<- #
        d_rec, d_disp, d_shift, pos_disp, pos_shift = _dd_

        print 'rec density:', d_rec.min(), d_rec.max(), d_rec.mean()
        print 'shift density:', d_shift.min(), d_shift.max(), d_shift.mean()
        print 'displaced density:', d_disp.min(), d_rec.max(), d_rec.mean()

	dt_rec=d_rec
	print 'rec delta:', dt_rec.min(), dt_rec.max(), dt_rec.mean()



        if p.save_data:

	    try: fn_write
	    except: pass

	    #np.savez(fn_write, pos_displaced=np.rollaxis(disp, 0, len(disp.shape)) )
	    np.savez(fn_write, d_rec=d_rec, d_disp=d_disp, d_shift=d_shift, \
	             pos_disp=pos_disp, pos_shift=pos_shift )


    elif p.import_reced_data==True:
        #->> 
	f=np.load(fn_write)
        d_rec, d_disp, d_shift = f['d_rec'], f['d_disp'], f['d_shift']




    ''' ->> analysis the data <<- '''
    do_powerspectrum = True

    if do_powerspectrum==True:
        kspace='linear'  #'log_linear'
    
        k_rec, pk_rec=psxi.pk(dt_rec, boxsize=p.boxsize, kspace=kspace)
        k_ori, pk_ori=psxi.pk(delta, boxsize=p.boxsize, kspace=kspace)
        k_disp, pk_disp=psxi.pk(d_disp, boxsize=p.boxsize, kspace=kspace)
        k_shift, pk_shift=psxi.pk(d_shift, boxsize=p.boxsize, kspace=kspace)
    
        if True:
            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
        
            ax[0].loglog(k_ori, pk_ori, 'r--')
            ax[0].loglog(k_rec, pk_rec, 'k-')
            ax[0].loglog(k_disp, pk_disp, 'b-')
            ax[0].loglog(k_shift, pk_shift, 'y-')
    
            ax[1].plot(k_ori, pk_rec/pk_ori, 'k-') 
            ax[1].plot(k_ori, pk_disp/pk_ori, 'r-') 
            ax[1].plot(k_ori, pk_shift/pk_ori, 'b:') 
            ax[1].set_xscale("log")
    
            fig.savefig(root+'figure/ps_disp_shift_comp.png')
            pl.show()


    """


    
    # ->> The End <<- #
    p.finalize()
