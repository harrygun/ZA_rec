import sys
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
    'cal_rect_transfer_func':    True,
    'disp_transfunc_fname':   'rk.dat',
    'raw_disp_field_fname':       'a.dat',
    'stat_disp_field_fname':       'a.dat',
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
    'py_stat_model_PDF':   False,
    'py_stat_potential_model_PDF':   True,
    }


if __name__=='__main__':

    # ->> initialization <<- #
    sec='Rect'
    init_dict=myDict(param_dict)+myDict(prog_control)
    p=pc.prog_init(section=sec, **init_dict)
    p.smooth_R = p.smooth_r

    root=p.folder

    do_simulation_test=False
    do_displacement_test=True

    if do_simulation_test:
        #->> 
        fn_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/node0/100.000xv0.dat'
	fn_pid='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/node0/PID0.ic'
	#fn_pid='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/node0/0.000PID0.dat'

	print 'read fn_pid ', fn_pid

        #->> position <<- #
        #pos, vel=mio.read_cita_simulation(fn_part, p.ngrid)
	#del vel
	#->> pid <<- #
        pid=mio.read_cita_simulation_pid(fn_pid, p.ngrid, head_size=4)

        #->>  <<-#
        print len(pid), pid.min(), pid.max()

	ll=np.arange(len(pid))+1
        err=ll-pid

	print 'len of non-zeros: ', len(np.where(err!=0)[0])


    if do_displacement_test:
        #->> 
	fn_disp='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/rec_data/stat_disp_0_100.dat'
        nblock=6
        dd=rd.rblock(p.stat_disp_field_fname, p.ngrid**3*nblock, dtype='float').reshape(nblock,p.ngrid,p.ngrid,p.ngrid)
	print dd.shape

        disp, disp_lpt=dd[:3], dd[3:]

        if True:
            # ->>  <<- #
	    bd=10

            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=1, 100
                
            cb1=ax[0].imshow(disp[axis,:,:,nsl])
            cb1=ax[1].imshow(disp_lpt[axis,:,:,nsl])

	    #pl.colorbar(cb1)
	    #pl.colorbar(cb2)

                
            pl.tight_layout()
            pl.show()





    # ->> The End <<- #
    p.finalize()

