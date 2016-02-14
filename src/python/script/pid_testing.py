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

    do_simulation_test=True

    if do_simulation_test:
        #->> 
        fn_part='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/node0/100.000xv0.dat'
	fn_pid='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml_pid/node0/PID0.ic'

	print 'read fn_pid ', fn_pid

        #->> position <<- #
        pos, vel=mio.read_cita_simulation(fn_part, p.ngrid)
	del vel
	#->> pid <<- #
        pid=mio.read_cita_simulation_pid(fn_pid, p.ngrid)

        #->>  <<-#
        print len(pid), pid.min(), pid.max()





    # ->> The End <<- #
    p.finalize()

