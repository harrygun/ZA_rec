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
    }


if __name__=='__main__':

    # ->> initialization <<- #
    sec='Rect'
    init_dict=myDict(param_dict)+myDict(prog_control)
    p=pc.prog_init(section=sec, **init_dict)
    p.smooth_R = p.smooth_r

    root=p.folder


    # ->> import data <<- #
    if (p.py_import_density_field==True):

        if p.import_format=='cita_simulation':
            p.rec_fname=p.reconstructed_fname+'_'+p.smooth_type+'_R'+str(p.smooth_R)+'.dat'
            print 'reading data ... ', p.rec_fname

            f_rec=rd.rblock(p.rec_fname, p.ngrid**3*3, dtype='float').reshape(3,p.ngrid,p.ngrid,p.ngrid)
            drec, d_disp, d_shift=f_rec

            d_ori=rd.rblock(p.original_density_fname, p.ngrid**3, dtype='float').reshape(p.ngrid,p.ngrid,p.ngrid)

            print 'density shape:', drec.shape, d_disp.shape, d_shift.shape, d_ori.shape
	    print 'density min/max:', drec.min(), drec.max(), d_disp.min(), d_disp.max(), d_shift.min(), d_shift.max()

        else:
            raise Exception



    # ->> power spectrum measurement <<- #
    if (p.py_pk_comparison==True):

        ki_ori, pki_ori=ps.pk(d_ori, boxsize=p.boxsize)

        ki_rec, pki_rec=ps.pk(drec, boxsize=p.boxsize)
        ki_disp, pki_disp=ps.pk(d_disp, boxsize=p.boxsize)
        ki_shift, pki_shift=ps.pk(d_shift, boxsize=p.boxsize)

        #ki_disp_ivf, pki_disp_ivf=ps.pk(d_disp_ivf, boxsize=p.boxsize)

        #- >>
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
        ax[0].loglog(ki_ori, pki_ori, 'y-')

	ax[1].semilogx(ki_ori, pki_shift/pki_ori, 'b:')

	pl.tight_layout()
	pl.show()






    # ->> The End <<- #
    p.finalize()

