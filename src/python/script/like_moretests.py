import matplotlib
import matplotlib.pyplot as pl
import matplotlib.colors as colors

import numpy as np
import array as arr
import pynbody as pn

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.myplot as mpl
import genscript.read as rd

import rect.misc.io as mio
import rect.misc.cyth.cgalio as cg
import rect.misc.cic as mcic
import rect.misc.ps as psor
import rect.fourier.psxi as ps



import rect.misc.file_import as fimp
import rect.misc.power_update as pu




"""
def disp_likelihood_test(p):

    print 'likelihood fname:', p.likelihood_test_fname
    z_init = 100.
    
    # ->> import testing data <<- #
    dd=rd.rblock(p.likelihood_test_fname, p.ngrid**3*7, \
                 dtype='float').reshape(7,p.ngrid,p.ngrid,p.ngrid)
    disp, disp_model = dd[:3,...], dd[3:,...]
    d=dd[-1,...]
    
    
    _disp_correlation_=True
    _disp_transfer_ = False
    
    if _disp_transfer_:
        # ->> calculate the transfer function of displacement field <<- #
        nplt, ncol = 3, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)
    
        lw1=['k-', 'r-', 'b-', 'g-']
        lw2=['k--', 'r--', 'b--', 'g--']
    
        for i in range(3):
            k1, pk1=psor.cross(disp[i], disp_model[i], boxsize=p.boxsize)
            k2, pk2=psor.pk(disp_model[i], boxsize=p.boxsize)
    
    	ax[i].plot(k1, pk2/pk1, lw1[i])
    
        pl.show()
    
    
    if _disp_correlation_:
        #->> calculate the correlation function between displacement <<- #
    
        nplt, ncol = 3, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)
    
        lw1=['k-', 'r-', 'b-', 'g-']
        lw2=['k--', 'r--', 'b--', 'g--']
    
        for i in range(3):
            k_, pk_=ps.pk(disp[i], boxsize=p.boxsize)
    	ax[i].loglog(k_, pk_, lw1[i])
            k_, pk_=ps.pk(disp_model[i], boxsize=p.boxsize)
    	ax[i].loglog(k_, pk_, lw2[i])
    
        pl.show()
    
    
    # ->> make plots <<- #
    if True:
        nplt, ncol = 3, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)
        n_bin=500
        color=['g', 'r', 'b', 'y']
    
        drange=[-20,20]
    
        for i in range(3):
            ax[i].hist(disp[i].flatten(), bins=n_bin, range=drange, \
                       histtype='step', color=color[0])
            ax[i].hist(disp_model[i].flatten(), bins=n_bin, \
                     range=drange, histtype='step', color=color[1])
            ax[i].hist((disp[i]+disp_model[i]).flatten(), bins=n_bin, \
                     range=drange, histtype='step', color=color[2])
            ax[i].hist((disp[i]-disp_model[i]).flatten(), bins=n_bin, \
                     range=drange, histtype='step', color=color[3])
    
        pl.tight_layout()
        pl.show()
    
    if True:
        # ->> show images <<- #
        nplt, ncol = 3, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
        	                      gap_size=0.5,return_figure=True)
        axis, nsl=0, 100
            
            #cb1=ax[0].imshow(disp[axis,:,:,nsl])
            cb1=ax[0].imshow(disp[axis,15:-15,15:-15,nsl])
            cb2=ax[1].imshow(-disp_model[axis,15:-15,15:-15,nsl])

            cb2=ax[2].imshow(d[15:-15,15:-15,nsl])

	    #fig.colorbar(cb1) #, orientation='horizontal')
	    #fig.colorbar(cb2) #, orientation='horizontal')

            
            pl.tight_layout()
            pl.show()




    return
"""





param_dict={
    'smooth_R_list_type':  'linear', 
    'smooth_R':    10,
    'smooth_r':    10,
    'smooth_type': 'Gaussian', 
    'boxsize':     512.,
    'ngrid':       256, 
    'redshift':    0.,
    'folder':      '~/',
    'import_format':   'cita_simulation',
    'original_density_fname':  'x.dat',
    'reconstructed_fname': 'y.dat',
    'other_test_fname':    'z.dat', 
    'particle_file_name':  'z.dat',
    'disp_transfunc_fname':   'rk.dat',
    'py_import_density_field':   True,
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    }

prog_control={
    #-------------------------------#
    'py_do_testing': False, 
    #-------------------------------#
    'py_image_comparison':  True,
    'py_pk_comparison':       False,
    'py_cf_comparison':       True,
    'py_part_position_check':  False,
    #-------------------------------#
    'do_likelihood_testing':   True,
    'likelihood_test_fname':   'x.dat',
    }


if __name__=='__main__':

    # ->> initialization <<- #
    sec='Rect'
    init_dict=myDict(param_dict)+myDict(prog_control)
    p=pc.prog_init(section=sec, **init_dict)
    p.smooth_R = p.smooth_r

    root=p.folder

    

    ''' ->> import data <<- '''
    print 'likelihood fname:', p.likelihood_test_fname
    z_init = 100.

    # ->> import testing data <<- #
    dd=rd.rblock(p.likelihood_test_fname, p.ngrid**3*7, \
	             dtype='float').reshape(7,p.ngrid,p.ngrid,p.ngrid)
    #disp, disp_model = dd[:3,...], dd[3:,...]
    bd=10
    disp, disp_model = dd[:3,bd:-bd,bd:-bd,bd:-bd], dd[3:,bd:-bd,bd:-bd,bd:-bd],

    di=dd[-1,...]

    #
    print 'mean:', disp.min(), disp_model.min()

    #->> final density <<- #
    df=rd.rblock(p.original_density_fname, p.ngrid**3, \
                 dtype='float').reshape(p.ngrid, p.ngrid, p.ngrid)


    # ->> define some other useful variables <<- #
    lw1=['k-', 'r-', 'b-', 'g-']


    # ->> get the divergence of the displacement field <<- #
    
    #->> 
    if False:
        print 'disp rk fname:', p.disp_transfunc_fname

        nplt, ncol = 1, 1
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
	                          gap_size=0.5,return_figure=True)

        lw1=['k-', 'r-', 'b-', 'g-']
        lw2=['k--', 'r--', 'b--', 'g--']
        
	cd_k, cd_p=[], []
        for i in range(3):

            k1, pk1=psor.cross(disp[i], disp_model[i], boxsize=p.boxsize)
            k2, pk2=psor.pk(disp_model[i], boxsize=p.boxsize)
            k3, pk3=psor.pk(disp[i], boxsize=p.boxsize)

	    cr=pk1/np.sqrt(pk3*pk2)

	    cd_k.append(k1)
	    cd_p.append(cr)

            ax[0].semilogx(k1, pk1/np.sqrt(pk3*pk2), lw1[i])
	    ax[0].set_ylim([0, 1.05])


        # ->> save some <<- #
	ddd=np.array([cd_k,cd_p])
	print 'R(k) shape:', ddd.shape
	ddd.tofile(p.disp_transfunc_fname)

        pl.show()



    if True:
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)

        nsl=100
	axis=0

        cb1=ax[0].imshow(disp[axis,:,:,nsl])
        cb2=ax[1].imshow(disp_model[axis,:,:,nsl])

	pl.colorbar(cb1)
	pl.colorbar(cb2)

        pl.show()




    # ->> The End <<- #
    p.finalize()

