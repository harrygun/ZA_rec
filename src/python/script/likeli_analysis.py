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
#import rect.misc.ps as ps
import rect.fourier.psxi as ps



import rect.misc.file_import as fimp
import rect.misc.power_update as pu




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





    if (p.do_likelihood_testing==True):
        print 'likelihood fname:', p.likelihood_test_fname

        # ->> import testing data <<- #
        dd=rd.rblock(p.likelihood_test_fname, p.ngrid**3*6, \
	             dtype='float').reshape(6,p.ngrid,p.ngrid,p.ngrid)
	disp, disp_model = dd[:3,...], dd[3:,...]

        # ->> make plots <<- #
	if True:
            nplt, ncol = 4, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
	                          gap_size=0.5,return_figure=True)
            n_bin=500
	    color=['g', 'r', 'b']

            for i in range(3):
                ax[0].hist(disp[i].flatten(), bins=n_bin, range=[-20,20], \
	                   histtype='step', color=color[i])
                ax[1].hist(disp_model[i].flatten(), bins=n_bin, \
	                 range=[-20,20], histtype='step', color=color[i])
                ax[2].hist((disp[i]+disp_model[i]).flatten(), bins=n_bin, \
	                 range=[-20,20], histtype='step', color=color[i])
                ax[3].hist((disp[i]-disp_model[i]).flatten(), bins=n_bin, \
	                 range=[-20,20], histtype='step', color=color[i])

	    pl.tight_layout()
	    pl.show()

        if True:
            # ->> show images <<- #
            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=1, 50
            
            #cb1=ax[0].imshow(disp[axis,:,:,nsl])
            cb1=ax[0].imshow(disp[axis,nsl,15:-15,15:-15])
            cb2=ax[1].imshow(disp_model[axis,15:-15,15:-15,nsl])

	    #fig.colorbar(cb1) #, orientation='horizontal')
	    #fig.colorbar(cb2) #, orientation='horizontal')

            
            pl.tight_layout()
            pl.show()







    # ->> The End <<- #
    p.finalize()
