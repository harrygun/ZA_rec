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

    

    ''' ->> import data <<- '''
    print 'likelihood fname:', p.likelihood_test_fname
    z_init = 100.

    # ->> import testing data <<- #
    dd=rd.rblock(p.likelihood_test_fname, p.ngrid**3*7, \
	             dtype='float').reshape(7,p.ngrid,p.ngrid,p.ngrid)
    disp, disp_model = dd[:3,...], dd[3:,...]
    di=dd[-1,...]

    #->> final density <<- #
    df=rd.rblock(p.original_density_fname, p.ngrid**3, \
                 dtype='float').reshape(p.ngrid, p.ngrid, p.ngrid)

    # define controllers:  
    _density_ZA_cor_ = False
    _density_propagator_ = True
    _density_map_cmp_ = False

    # ->> define some other useful variables <<- #
    lw1=['k-', 'r-', 'b-', 'g-']


    # ->> get the divergence of the displacement field <<- #

    if _density_propagator_:
        #->> 
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
	                          gap_size=0.5,return_figure=True)

	di=di[1:,1:,1:]
	df=df[1:,1:,1:]

        k1, ck1=psor.cross(di, df, boxsize=p.boxsize)
        k2, pk2=psor.pk(di, boxsize=p.boxsize)
        k3, pk3=psor.pk(df, boxsize=p.boxsize)

	#print 'len:', k1.shape, ck1.shape, k2.shape, pk2.shape, k3.shape, pk3.shape

	ax[0].semilogx(k1, ck1/np.sqrt(pk2[1:]*pk3[1:]))
	#ax[0].semilogx(k1, ck1/np.sqrt(pk2*pk3))

	ax[1].loglog(k1, ck1,'k-')
	ax[1].loglog(k2, pk2,'r-')
	ax[1].loglog(k3, pk3,'b-')

	pl.show()


    if _density_map_cmp_:
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)


        nsl=20
        #cb=ax[0].imshow(di[:,:,nsl])

	print 'di:', di[0,0,0]

        ax[0].imshow(di[1:,1:,nsl])
        ax[1].imshow(df[:,:,nsl]-df.min()+1e-3, norm=colors.LogNorm())


        pl.show()


    if _density_ZA_cor_:
        nplt, ncol = 1, 1
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
	                          gap_size=0.5,return_figure=True)
        tt=np.zeros(di.shape)

        for i in range(3):
	  #tt+=disp[i]*disp[i]
	  tt+=disp_model[i]*disp_model[i]

	tt=np.sqrt(tt)

        k1, ck1=psor.cross(tt, di, boxsize=p.boxsize)
        k2, pk2=psor.pk(di, boxsize=p.boxsize)


        ax[0].semilogx(k1, ck1/pk2, lw1[i])
        pl.tight_layout()
        pl.show()






    # ->> The End <<- #
    p.finalize()

