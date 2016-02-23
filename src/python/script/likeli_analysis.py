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
    'cal_rect_transfer_func':    False,
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
    'py_stat_potential_model_PDF':   False,
    }


if __name__=='__main__':

    # ->> initialization <<- #
    sec='Rect'
    init_dict=myDict(param_dict)+myDict(prog_control)
    p=pc.prog_init(section=sec, **init_dict)
    p.smooth_R = p.smooth_r

    root=p.folder



    if (p.cal_rect_transfer_func==True):
        nblock=6
        dd=rd.rblock(p.raw_disp_field_fname, p.ngrid**3*nblock, dtype='float').reshape(nblock,p.ngrid,p.ngrid,p.ngrid)

	#->> discard boundary data <<- #
        bd=5
	bsize=p.boxsize-2.*bd*p.boxsize/float(p.ngrid)
	#bsize=p.boxsize
	print 'bsize:', bsize

        disp, disp_model = dd[:3,bd:-bd,bd:-bd,bd:-bd], dd[3:,bd:-bd,bd:-bd,bd:-bd],
        #disp, disp_model = dd[:3], dd[3:]

        if True:
            nplt, ncol = 2, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=1, 20
                
            cb1=ax[0].imshow(disp[axis,:,:,nsl])
            cb2=ax[1].imshow(disp_model[axis,:,:,nsl])


	    #pl.colorbar(cb1)
	    pl.colorbar(cb2)

                
            pl.tight_layout()
            pl.show()


	quit()



        lw1=['k-', 'r-', 'b-', 'g-']
        lw2=['k--', 'r--', 'b--', 'g--']

        nplt, ncol = 1, 1
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
	                          gap_size=0.5,return_figure=True)

	_cd_k, _cd_p=[], []
        for i in range(3):
            k1, pk1=psor.cross(disp[i], disp_model[i], boxsize=bsize)
            k2, pk2=psor.pk(disp_model[i], boxsize=bsize)
            k3, pk3=psor.pk(disp[i], boxsize=bsize)

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

            ax[0].semilogx(k1, pk1/np.sqrt(pk3*pk2), lw1[j])
	    ax[0].set_ylim([0, 1.05])

        # ->> close <<- #
	f.close()

        pl.show()




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
        '''
        nb1, nb2=9, 10
	ntrim=10

	ng=p.ngrid
	tng=p.ngrid-2*ntrim

        dd=rd.rblock(p.stat_disp_field_fname, ng**3*nb1+tng**3*nb2, dtype='float')
        dd1=dd[:ng**3*nb1].reshape(nb1,ng,ng,ng)
        dd2=dd[ng**3*nb1:].reshape(nb2,tng,tng,tng)

        bd=ntrim
        #disp = dd1[:3,bd:-bd,bd:-bd,bd:-bd]
        #disp_lpt = dd1[3:6,bd:-bd,bd:-bd,bd:-bd]
	#disp_mc=dd1[6:,bd:-bd,bd:-bd,bd:-bd]

        disp = dd1[:3]
        disp_lpt = dd1[3:6]
	disp_mc=dd1[6:]

        div, phi = dd2[0], dd2[1]
        disp_phi=dd2[2:5]
        div_lpt,phi_lpt =dd2[5], dd2[6]
	disp_phi_lpt=dd2[7:]
        '''

        nb=19
	ntrim=5
	tng=p.ngrid-2*ntrim
        dd=rd.rblock(p.stat_disp_field_fname, tng**3*nb, dtype='float').reshape(nb,tng,tng,tng)

        disp, disp_lpt, disp_mc = dd[:3], dd[3:6], dd[6:9]

        div, phi = dd[9+0], dd[9+1]
        disp_phi=dd[9+2:9+5]
        div_lpt,phi_lpt =dd[9+5], dd[9+6]
	disp_phi_lpt=dd[9+7:]

        ''' ->> end of data importing <<- '''

	disp_phimc=disp_phi-disp_lpt
	disp_rot=disp-disp_phi

	#->>
	print disp.shape, div.shape, phi.shape, disp_phi.shape, 
	print disp_lpt.shape, div_lpt.shape, phi_lpt.shape, disp_mc.shape


	err=(disp_phi_lpt-disp_lpt)
	print 'err:', err.min(), err.max()


	if False:
	    # ->> 2D histogram <<- #
            nplt, ncol = 6, 3
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                                  gap_size=0.5,return_figure=True)
            n_bin=200

            drange_div=[[-3,3], [-3,3]]
            drange_phi=[[-10,10], [-200,200]]

            drange_oo=[[-10,10], [-15,15]]

            ax[0].hist2d((div-div_lpt).flatten(), div_lpt.flatten(), 
                           bins=n_bin, range=drange_div, normed=True)

            ax[1].hist2d((phi-phi_lpt).flatten(), phi_lpt.flatten(), 
                           bins=n_bin, normed=True, range=drange_phi) 

            ax[2].hist2d(phi.flatten(), phi_lpt.flatten(), 
                           bins=n_bin, normed=True )#, range=drange_oo) 

            ax[3].hist2d(div_lpt.flatten(), div.flatten(), 
                           bins=n_bin, normed=True )#, range=drange_oo) 


            drange_rot=[[-1,1], [-200,200]]
            #for i in range(3):
            #    ax[i+3].hist2d(disp_lpt[i].flatten(), disp_phi.flatten(), 
            #               bins=n_bin, normed=True, range=drange_rot)


            pl.tight_layout()
            pl.show()




        if False:
	    # ->> 1D histogram <<- #
            nplt, ncol = 3, 3
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                                  gap_size=0.5,return_figure=True)
            n_bin=500
            color=['g', 'r', 'b', 'y', 'k', 'm']
    
            #drange=[-10,10]
            drange=[-80,80]
            #drange=[-20,20]
    
            for i in range(3):
                ax[i].hist(disp[i].flatten(), bins=n_bin, range=drange, \
                         normed=True, histtype='step', color=color[0])
                ax[i].hist(disp_lpt[i].flatten(), bins=n_bin, \
                         normed=True, range=drange, histtype='step', color=color[1])
                ax[i].hist((disp_mc[i]).flatten(), bins=n_bin, \
                        normed=True, range=drange, histtype='step', color=color[2])
                ax[i].hist((disp_phi[i]).flatten(), bins=n_bin, \
                        normed=True, range=drange, histtype='step', color=color[3])
                ax[i].hist((disp_phimc[i]).flatten(), bins=n_bin, \
                        normed=True, range=drange, histtype='step', color=color[4])
                #ax[i].hist((disp_rot[i]).flatten(), bins=n_bin, \
                #        normed=True, range=drange, histtype='step', color=color[5])

                ax[i].hist((err[i]).flatten(), bins=n_bin, \
                        normed=True, range=drange, histtype='step', color=color[5])

            pl.tight_layout()
            pl.show()



        if True:
            # ->>  <<- #

            nplt, ncol = 4, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
            	                      gap_size=0.5,return_figure=True)
            axis, nsl=0, 20
                
            cb1=ax[0].imshow(disp[axis,:,:,nsl])
            #cb2=ax[1].imshow(disp_lpt[axis,:,:,nsl])
            cb2=ax[1].imshow(disp_phi[axis,:,:,nsl])

            cb3=ax[2].imshow(disp_phi_lpt[axis,:,:,nsl])
            cb4=ax[3].imshow(disp_lpt[axis,:,:,nsl])


            #cb3=ax[2].imshow(phi[...,nsl])
            #cb4=ax[3].imshow(phi_lpt[...,nsl])

            #cb3=ax[2].imshow(div[...,nsl])
            #cb4=ax[3].imshow(div_lpt[...,nsl])

	    pl.colorbar(cb1)
	    pl.colorbar(cb2)
	    pl.colorbar(cb3)
	    pl.colorbar(cb4)

                
            pl.tight_layout()
            pl.show()



        if False:
	    # ->> test boundary <<- #
            nplt, ncol = 3, 3
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                                  gap_size=0.5,return_figure=True)
            n_bin=500
            color=['g', 'r', 'b', 'y', 'k', 'm']
    
            #drange=[-20,20]

            bd_list=[20, 15, 10, 5]
            #bd=10
            for i in range(3):
	        for j in range(len(bd_list)):
		    bd=bd_list[j]
                    d_ = np.roll(np.roll(np.roll(disp[i],bd,axis=0),bd,axis=1),bd,axis=2)[0:2*bd,0:2*bd,0:2*bd].flatten()
                    ax[i].hist(d_, bins=n_bin, histtype='step', color=color[j]) #range=drange, \

            pl.tight_layout()
            pl.show()





    # ->> The End <<- #
    p.finalize()

