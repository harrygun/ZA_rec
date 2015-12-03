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
import rect.misc.ps as ps


import rect.misc.file_import as fimp






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'smooth_R': 15.,
    'smooth_type': 'Gaussian', 
    'smooth_R_list_type':  'linear', 
    'boxsize': 512.,
    'ngrid':    256, 
    'redshift':   0.,
    'folder':     '~/',
    'import_format':   'cita_simulation',
    'original_density_fname':  'x.dat',
    'reconstructed_fname': 'y.dat',
    }

prog_control={
    #-------------------------------#
    'py_do_testing': False, 
    #-------------------------------#
    'py_image_comparison':  True,
    'py_pk_comparison':       False,
    #-------------------------------#
    }


if __name__=='__main__':

    # ->> initialization <<- #
    sec='Rect'
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(section=sec, **init_dict)

    root=p.folder


    # ->> import data <<- #
    if p.import_format=='cita_simulation':

        print 'reading data ... '

        f_rec=rd.rblock(p.reconstructed_fname, p.ngrid**3*3, dtype='float').reshape(3,p.ngrid,p.ngrid,p.ngrid)
	drec, d_disp, d_shift=f_rec

	d_ori=rd.rblock(p.original_density_fname, p.ngrid**3, dtype='float').reshape(p.ngrid,p.ngrid,p.ngrid)

	print 'density shape:', drec.shape, d_disp.shape, d_shift.shape, d_ori.shape

    else:
        raise Exception

    


    # ->> power spectrum measurement <<- #
    if (p.py_pk_comparison==True):

        ki_ori, pki_ori=ps.pk(d_ori, boxsize=p.boxsize)

        ki_rec, pki_rec=ps.pk(drec, boxsize=p.boxsize)
        ki_disp, pki_disp=ps.pk(d_disp, boxsize=p.boxsize)
        ki_shift, pki_shift=ps.pk(d_shift, boxsize=p.boxsize)

        #- >>
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
        ax[0].loglog(ki_ori, pki_ori, 'y-')
        ax[0].loglog(ki_rec, pki_rec, 'k-')
        ax[0].loglog(ki_disp, pki_disp, 'r--')
        ax[0].loglog(ki_shift, pki_shift, 'b:')

	ax[1].semilogx(ki_ori, pki_rec/pki_ori, 'k-')
	ax[1].semilogx(ki_ori, pki_disp/pki_ori, 'r--')
	ax[1].semilogx(ki_ori, pki_shift/pki_ori, 'b:')

	pl.show()



    #->> image  comparison <<- #
    if (p.py_image_comparison==True):
        print 'making plots ... '

        nplt, ncol = 3, 3
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)

	sl=100
        cb0=ax[0].imshow(drec[:,:,sl]-drec.min()+1e-3, norm=colors.LogNorm())
        cb1=ax[1].imshow(d_disp[:,:,sl]-d_disp.min()+1e-3, norm=colors.LogNorm())
        cb2=ax[2].imshow(d_shift[:,:,sl]-d_shift.min()+1e-3, norm=colors.LogNorm())

	#pl.colorbar(cb0)
	#pl.colorbar(cb1)
	#pl.colorbar(cb2)

        #fig.savefig(root+'figure/ps_comp.png')
	pl.show()




    # ->> The End <<- #
    p.finalize()

