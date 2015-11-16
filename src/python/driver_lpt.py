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





def Testing_portal(p, data, import_data_type='field'):

    if False:
        if import_data_type=='field':
            # ->> 
	    dmean=np.mean(data)
	    delta = data/dmean -1.
            #phi, phi_ij=ptt.Poisson3d(delta, boxsize=p.boxsize, return_hessian=True, \
            #                          smooth_R=None, smooth_type=None)
            phi, phi_ij=ptt.Poisson3d(delta, boxsize=p.boxsize, return_hessian=True, \
                                      smooth_R=p.smooth_R, smooth_type=p.smooth_type)
            # ->> 
            d=dmean*(1+phi_ij[0,0]+phi_ij[1,1]+phi_ij[2,2])
    	    print 'd shape:', d.shape

	    # ->> 
	    da = dmean*(1.+ptt.Laplacian(phi, boxsize=p.boxsize, Lap_type=1))
	    print 'da shape:', da.shape
    
            if True:
                nplt = 3
                ncol = 3
                fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=2.,gap_size=0.15,return_figure=True)

    	        dd1, dd2, dd3=data[100,:,:], d[100,:,:], da[100,:,:]
    	        print 'min, max:', dd1.min(), dd1.max(), dd2.min(), dd2.max(), dd3.min(), dd3.max()
    
                ax[0].imshow(np.flipud(dd1), norm=colors.LogNorm(vmin=dd1.min(), vmax=dd1.max()) )
                ax[1].imshow(np.flipud(dd2), norm=colors.LogNorm(vmin=dd2.min(), vmax=dd2.max()) )
                ax[2].imshow(np.flipud(dd3), norm=colors.LogNorm(vmin=dd3.min(), vmax=dd3.max()) )
    
    	    pl.show()


    if True:

        if import_data_type=='field':
	    dmean=np.mean(data)
	    delta = data/dmean -1.

            #phi=ptt.Poisson3d(delta, boxsize=p.boxsize, smooth_R=None, smooth_type=None)
            phi=ptt.Poisson3d(delta, boxsize=p.boxsize, smooth_R=p.smooth_R, smooth_type=p.smooth_type)

	    if True:
	        # ->> testing <<- #
                gd=np.array(np.gradient(phi))
	        print gd.shape

                if True:
                    nplt = 6
                    ncol = 3
                    fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=2.,gap_size=0.15,return_figure=True)
                    ax[0].imshow(np.flipud(phi[100,:,:]))

	            for i in range(3):
		        dd_ = np.fabs(gd[i,100,:,:])
                        ax[i+1].imshow(np.flipud(dd_), norm=colors.LogNorm(vmin=dd_.min(), vmax=dd_.max()) )

                    _dd = gd[0,100,:,:]**2.+gd[1,100,:,:]**2.+gd[2,100,:,:]**2.
                    ax[nplt-2].imshow(np.flipud(_dd), norm=colors.LogNorm(vmin=_dd.min(), vmax=_dd.max()) )

                    dd1=data[100,:,:]
                    ax[nplt-1].imshow(np.flipud(dd1), norm=colors.LogNorm(vmin=dd1.min(), vmax=dd1.max()) )

	            pl.show()


    if False:
        # ->> testing <<- #
        pl.plot(y[:,:,100]*p.nbin, z[:,:,100]*p.nbin, 'k.', alpha=0.3)
    
        #_dd=np.swapaxes(den[100,:,:], 0, 1)
        #_dd=den[100,:,:]
        _dd=np.flipud(den[100,:,:])
        #pl.pcolormesh(_dd, norm=colors.LogNorm(vmin=_dd.min(), vmax=_dd.max()), alpha=0.7)
        pl.imshow(_dd, norm=colors.LogNorm(vmin=_dd.min(), vmax=_dd.max()), alpha=0.8)
    
        pl.show()
    return










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
    import_type='all'
    #import_type='field'

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

        dd=fimp.import_cita_simulation(p, fn_part, fn_field, import_data_type=import_type)

	# ->> estimate mass resolution <<- #

        p.particle_mass = mcic.mass_resolution(p, z=0., boxsize_unit='Mpc/h')
	print 'particle mass:', p.particle_mass

    if p.import_format=='cita_simulation_highres':
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



    ''' -------------------------------------------------
               ->>      do some testing      <<- 
	-------------------------------------------------     
    '''
    if p.do_testing==True:
        Testing_portal(p, dd, import_data_type=import_type)


    ''' ->> Performing Lagrangian reconstruction <<- '''
    if p.do_LPT_rec==True:
        pos, delta =dd
	print 'data shape:', pos.shape, delta.shape, delta.min(), delta.max(), delta.mean()



        ''' ->> perform Lagrangian Reconstruction <<-  '''

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
    
        k_rec, pk_rec=ps.pk(dt_rec, boxsize=p.boxsize)
        k_ori, pk_ori=ps.pk(delta, boxsize=p.boxsize)

        k_disp, pk_disp=ps.pk(d_disp, boxsize=p.boxsize)
    
        if True:
            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
        
            ax[0].loglog(k_ori, pk_ori, 'r--')
            ax[0].loglog(k_rec, pk_rec, 'k-')
            ax[0].loglog(k_disp, pk_disp, 'b-')
    
            ax[1].plot(k_ori, pk_rec/pk_ori, 'k-') 
            ax[1].plot(k_ori, pk_disp/pk_ori, 'r-') 
            ax[1].set_xscale("log")
    
            fig.savefig(root+'figure/ps_disp_shift_comp.png')
            pl.show()


	# ->> correlation function <<- #
        #r_ori, xi_ori= 
    
    
    
    if True:
        # ->> comparison plot <<- #
    
        nplt = 2
        ncol = 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=10.,gap_size=0.5,return_figure=True)
    
        #ax[0].imshow(np.flipud(phi[100,:,:]))
        #ax[0].plot(pos[1,:,:,100], pos[2,:,:,100], 'k.', alpha=0.3)
        #ax[1].plot(disp[1,:,:,100], disp[2,:,:,100], 'r.', alpha=0.3)
    
        dat_rec=dt_rec[:,:,100]-1.01*np.min(dt_rec[:,:,100])
        dat_ori=delta[:,:,100]-1.01*np.min(delta[:,:,100])
        ax[0].imshow(np.flipud(dat_rec), norm=colors.LogNorm(vmin=dat_rec.min(),vmax=dat_rec.max()))
        ax[1].imshow(np.flipud(dat_ori), norm=colors.LogNorm(vmin=dat_ori.min(),vmax=dat_ori.max()))
    
        pl.show()
        #fig.savefig('rect.png')
    
    
    
    # ->> if making plots <<- #
    if False:
        nplt = 4
        ncol = 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=2.,gap_size=0.15,return_figure=True)
        ax[0].imshow(np.flipud(phi[100,:,:]))
    
        for i in range(3):
            ax[i+1].imshow(np.flipud(si[i,100,:,:]))
    
    
        pl.show()






    
    # ->> The End <<- #
    p.finalize()
