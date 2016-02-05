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

#import likeli_test as ltst







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



    if (p.cal_rect_transfer_func):
        nblock=7
        dd=rd.rblock(p.raw_disp_field_fname, p.ngrid**3*nblock, dtype='float').reshape(nblock,p.ngrid,p.ngrid,p.ngrid)

	#->> discard boundary data <<- #
        bd=10
        disp, disp_model = dd[:3,bd:-bd,bd:-bd,bd:-bd], dd[3:,bd:-bd,bd:-bd,bd:-bd],


	_cd_k, _cd_p=[], []
        for i in range(3):
            k1, pk1=psor.cross(disp[i], disp_model[i], boxsize=p.boxsize)
            k2, pk2=psor.pk(disp_model[i], boxsize=p.boxsize)
            k3, pk3=psor.pk(disp[i], boxsize=p.boxsize)

	    cr=pk1/np.sqrt(pk3*pk2)

	    _cd_k.append(k1)
	    _cd_p.append(cr)
	

	cd_k=np.array(_cd_k)
	cd_p=np.array(_cd_p)

        # ->> save data <<- #
        f=open(p.disp_transfunc_fname, "w")

	for i in range(cd_k.shape[-1]):
	    for j in range(3):
	        strr="{0}  {1}  ".format(cd_k[j,i],cd_p[j,i]).rstrip('\n')
                f.write(strr)
	    f.write("\n")


	f.close()
        # ->> close <<- #



    if (p.do_likelihood_testing==True):
        pass





    # ->> The End <<- #
    p.finalize()

