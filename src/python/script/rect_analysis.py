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
            print 'reading data ... '
            p.rec_fname=p.reconstructed_fname+'_'+p.smooth_type+'_R'+str(p.smooth_R)+'.dat'
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

	pl.tight_layout()
	pl.show()



    # ->> correlation function comparison <<- //
    if (p.py_cf_comparison==True):

        r_ori, xi_ori=pu.corfunk(d_ori, boxsize=p.boxsize, binsize=1)

        r_rec,   xi_rec=pu.corfunk(drec, boxsize=p.boxsize)
        r_disp,  xi_disp=pu.corfunk(d_disp, boxsize=p.boxsize)
        r_shift, xi_shift=pu.corfunk(d_shift, boxsize=p.boxsize)


        #- >>
        nplt, ncol = 1, 1
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.5,return_figure=True)
        ax[0].loglog(r_ori,  r_ori**2.*xi_ori, 'y-')
        ax[0].loglog(r_rec,  r_rec**2.*xi_rec, 'k-')
        ax[0].loglog(r_disp, r_disp**2.*xi_disp, 'r--')
        ax[0].loglog(r_shift,r_shift**2.*xi_shift, 'b:')

	#ax[1].semilogx(ki_ori, xi_rec/xi_ori, 'k-')
	#ax[1].semilogx(ki_ori, xi_disp/xi_ori, 'r--')
	#ax[1].semilogx(ki_ori, xi_shift/xi_ori, 'b:')

	pl.tight_layout()
	pl.show()




    #->> image  comparison <<- #
    if (p.py_image_comparison==True):
        print 'making plots ... '

        nplt, ncol = 4, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,gap_size=0.1,return_figure=True)

	sl=100
        cb0=ax[0].imshow(drec[:,:,sl]-drec.min()+1e-3, norm=colors.LogNorm())
        cb1=ax[1].imshow(d_disp[:,:,sl]-d_disp.min()+1e-3, norm=colors.LogNorm())
        cb2=ax[2].imshow(d_shift[:,:,sl]-d_shift.min()+1e-3, norm=colors.LogNorm())
        cb3=ax[3].imshow(d_ori[:,:,sl]-d_ori.min()+1e-3, norm=colors.LogNorm())


        #fig.savefig(root+'figure/ps_comp.png')
	pl.tight_layout()
	pl.show()


    # ->> check particles <<- #
    if (p.py_part_position_check==True):
        print 'checking moved particle position...'

        if False:
	    #->> check boundary of original particles <<- #
            pos, vel=mio.read_cita_simulation(p.particle_file_name, p.ngrid)
	    del vel
            for i in range(3):
                print 'pos:', pos[...,i].min(), pos[...,i].max()
	    quit()


	#->> import particles <<- #
        f=rd.rblock(p.other_test_fname, p.ngrid**3*3, dtype='float').reshape(2,p.ngrid**3*3)
        pos=f[0].reshape(p.ngrid,p.ngrid,p.ngrid,3)
        si=f[1].reshape(3,p.ngrid,p.ngrid,p.ngrid)

        for i in range(3):
	    print 'pos:', i, pos[...,i].min(), pos[...,i].max()
	    print 'si:', i, si[i,...].min(), si[i,...].max()

        nplt, ncol = 1, 1
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=8.,gap_size=0.1,return_figure=True)

	sl=100
	ax[0].plot(pos[:,:,sl,2], pos[:,:,sl,1], 'k.')

	pl.tight_layout()
	pl.show()




    # ->> The End <<- #
    p.finalize()

