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



    if (p.cal_rect_transfer_func==True):
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



    if (p.py_stat_model_PDF==True):

        nblock=9
        dd=rd.rblock(p.stat_disp_field_fname, p.ngrid**3*nblock, dtype='float').reshape(nblock,p.ngrid,p.ngrid,p.ngrid)

        #disp, disp_lpt, disp_mc=dd[:3,...], dd[3:6,...], dd[6:,...]
        bd=10
        disp = dd[:3,bd:-bd,bd:-bd,bd:-bd]
	disp_lpt=dd[3:6,bd:-bd,bd:-bd,bd:-bd]
	disp_mc=dd[6:,bd:-bd,bd:-bd,bd:-bd]

	print 'disp shape:', disp.shape, disp_lpt.shape, disp_mc.shape

	print disp.min(), disp.max()
	print disp_lpt.min(), disp_lpt.max()
	print disp_mc.min(), disp_mc.max()



        if True:
	    # ->> 1D histogram <<- #
            nplt, ncol = 3, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                                  gap_size=0.5,return_figure=True)
            n_bin=500
            color=['g', 'r', 'b', 'y']
    
            drange=[-10,10]
    
            for i in range(3):
                ax[i].hist(disp[i].flatten(), bins=n_bin, range=drange, \
                         normed=True, histtype='step', color=color[0])
                ax[i].hist(disp_lpt[i].flatten(), bins=n_bin, \
                         normed=True, range=drange, histtype='step', color=color[1])
                ax[i].hist((disp_mc[i]).flatten(), bins=n_bin, \
                        normed=True, range=drange, histtype='step', color=color[2])

            pl.tight_layout()
            pl.show()


        if True:
	    # ->>  2D histogram <<- #
            nplt, ncol = 3, 3
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                                  gap_size=0.5,return_figure=True)
            n_bin=500
            color=['g', 'r', 'b', 'y']
            drange=[[-3,3], [-8,8]]

            for i in range(3):
                ax[i].hist2d(disp_mc[i].flatten(), disp_lpt[i].flatten(), 
                           bins=n_bin, range=drange, normed=True)

            pl.tight_layout()
            pl.show()


        if False:
            # ->>  <<- #

            nplt, ncol = 3, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=0, 100
                
            cb1=ax[0].imshow(disp[axis,15:-15,15:-15,nsl])
            cb2=ax[1].imshow(disp_lpt[axis,15:-15,15:-15,nsl])
            cb3=ax[2].imshow(disp_mc[axis,15:-15,15:-15,nsl])

	    pl.colorbar(cb1)
	    pl.colorbar(cb2)
	    pl.colorbar(cb3)

                
            pl.tight_layout()
            pl.show()


    if (p.py_stat_potential_model_PDF==True):
        # ->> import data <<- #
        nblock=16
        dd=rd.rblock(p.stat_disp_field_fname, p.ngrid**3*nblock, dtype='float').reshape(nblock,p.ngrid,p.ngrid,p.ngrid)

        #disp, disp_lpt, disp_mc=dd[:3,...], dd[3:6,...], dd[6:,...]
        bd=10
        disp = dd[:3,bd:-bd,bd:-bd,bd:-bd]
        div  = dd[3,bd:-bd,bd:-bd,bd:-bd]
        phi  = dd[4,bd:-bd,bd:-bd,bd:-bd]
        disp_phi = dd[5:8,bd:-bd,bd:-bd,bd:-bd]

        disp_lpt = dd[8:11,bd:-bd,bd:-bd,bd:-bd]
        div_lpt  = dd[11,bd:-bd,bd:-bd,bd:-bd]
        phi_lpt  = dd[12,bd:-bd,bd:-bd,bd:-bd]

	disp_mc=dd[13:,bd:-bd,bd:-bd,bd:-bd]

	#->>
	#print disp.shape, div.shape, phi.shape, disp_phi.shape, 
	#print disp_lpt.shape, div_lpt.shape, phi_lpt.shape, disp_mc.shape

        #->> 

        print 'div min/max:', div.min(), div.max()
        print 'div_lpt min/max:', div_lpt.min(), div_lpt.max()

        print 'phi min/max:', phi.min(), phi.max()
        print 'phi_lpt min/max:', phi_lpt.min(), phi_lpt.max()

        if True:
            # ->>  <<- #

            nplt, ncol = 4, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=0, 20
                
            cb1=ax[0].imshow(disp[axis,15:-15,15:-15,nsl])
            cb2=ax[1].imshow(disp_lpt[axis,15:-15,15:-15,nsl])
            #cb3=ax[2].imshow(disp_phi[axis,15:-15,15:-15,nsl])

            cb3=ax[2].imshow(phi[15:-15,15:-15,nsl])
            cb4=ax[3].imshow(phi_lpt[15:-15,15:-15,nsl])

	    pl.colorbar(cb1)
	    pl.colorbar(cb2)
	    pl.colorbar(cb3)
	    pl.colorbar(cb4)

                
            pl.tight_layout()
            pl.show()



    # ->> The End <<- #
    p.finalize()

