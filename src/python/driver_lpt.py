import os
import numpy as np
import pylab as pl
#import pynbody as pn
import matplotlib.colors as colors

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.read as read
import genscript.myplot as mpl

import lagrangian.lpt_rec as lrec
import misc.io as mio
import fourier.potential as ptt





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

    return








'''-------------------------------------------------------------------------- '''

def import_data(p, droot_part, droot_field, fn, import_data_type='all'):

    print ' >> Importing data ...'

    if import_data_type=='all':
        datype=['particle', 'field']
    else:
        datype=[import_data_type]

    if 'particle' in datype:
        npart=int(p.nbin**3)
        _x, _y, _z = mio.read_cgal(droot_part, fn, npart, import_type='position')
        print _x.shape

        x=_x.reshape(p.nbin, p.nbin, p.nbin)
        y=_y.reshape(p.nbin, p.nbin, p.nbin)
        z=_z.reshape(p.nbin, p.nbin, p.nbin)
        pos=np.array([x, y, z])
    
    # ->> import field <<- #
    if 'field' in datype:
        fn_d = droot_field+fn+'.'+str(p.nbin)+'.fvol'
        den=read.rgrid(fn_d, ngrid=p.nbin, dtype='float', comp=1)
        print 'imported density shape:', den.shape

    if import_data_type=='all':
        return pos, den
    elif import_data_type=='particle':
        return pos
    elif import_data_type=='field':
        return den






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 1.,
    'smooth_type': 'Gaussian', 
    'smooth_R_list_type':  'linear', 
    'boxsize': 32.,
    'nbin':    256, 
    }

prog_control={
    #-------------------------------#
    'do_LPT_rec': False,
    #-------------------------------#
    'do_testing': True, 
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/result/'
    droot_part='/mnt/scratch-lustre/xwang/data/velinv/MIP/particle/'
    droot_field='/mnt/scratch-lustre/xwang/data/velinv/MIP/256/raw/'

    ''' -------------------------------------------------
	-------------------------------------------------     
    '''
    #import_type='all'
    import_type='field'

    # ->> imporot particles <<- #
    p.nbin=256
    fn='REAL_00-00-00'
    dd=import_data(p, droot_part, droot_field, fn, import_data_type=import_type)


    ''' -------------------------------------------------
               ->>      do some testing      <<- 
	-------------------------------------------------     
    '''
    if p.do_testing==True:
        Testing_portal(p, dd, import_data_type=import_type)


    ''' ->> Performing Lagrangian reconstruction <<- '''
    if p.do_LPT_rec==True:
        pos, den=dd
        dmean=np.mean(den)
        delta = den/dmean -1.
	print 'data shape:', pos.shape, den.shape, delta.shape

        if False:
	    # ->> testing <<- #
	    pl.plot(y[:,:,100]*p.nbin, z[:,:,100]*p.nbin, 'k.', alpha=0.3)

	    #_dd=np.swapaxes(den[100,:,:], 0, 1)
	    #_dd=den[100,:,:]
	    _dd=np.flipud(den[100,:,:])
            #pl.pcolormesh(_dd, norm=colors.LogNorm(vmin=_dd.min(), vmax=_dd.max()), alpha=0.7)
            pl.imshow(_dd, norm=colors.LogNorm(vmin=_dd.min(), vmax=_dd.max()), alpha=0.8)

	    pl.show()


        ''' ->> perform Lagrangian Reconstruction <<-  '''
        phi, si = lrec.lag_rec(p, pos, delta, boxsize=p.boxsize, smooth_R=p.smooth_R, \
                               smooth_type=p.smooth_type)

        # ->> if making plots <<- #
        if True:
            nplt = 4
            ncol = 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=2.,gap_size=0.15,return_figure=True)
            ax[0].imshow(np.flipud(phi[100,:,:]))

	    for i in range(3):
                ax[i+1].imshow(np.flipud(si[i,100,:,:]))


	    pl.show()







    
    # ->> The End <<- #
    p.finalize()
