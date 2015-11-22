import os
import numpy as np
#import pylab as pl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
#import pynbody as pn
import genscript.myplot as mpl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.interpolation as itp
import genscript.myutlis as mut
import genscript.myarray as mar

import fourier.potential as ptt
import misc.cic as mcic
import fourier.psxi as psxi
import misc.cyth.partmoving as pmv



displace_interpolation=True


def get_ZA_displacement(p, m, smooth_R=None, smooth_type=None):
    # ->> Poisson solver <<- #

    phi, phi_i=ptt.Poisson3d(m, boxsize=p.boxsize, smooth_R=smooth_R, \
                             smooth_type=smooth_type, return_gradient=True)

    if False:
        k_phi, pk_phi=psxi.pk(phi, boxsize=p.boxsize, kspace='linear')
        k_d, pk_d=psxi.pk(m, boxsize=p.boxsize, kspace='linear')

	pl.loglog(k_d, pk_d, 'r-')
	pl.loglog(k_phi, k_phi**4.*pk_phi, 'k-')
	pl.show()

        quit()

    if False:
        # ->> making plots <<- #
        pl.imshow(phi[:,:,100])
        pl.show()

    # ->> inverse of ZA: Psi_i = -partial_i * Laplacian^{-1} delta(tau) <<- #
    #dl=p.boxsize/np.float(p.nbin)
    #si = np.array(np.gradient(phi))/dl
    #print 'dl=', dl, '|si|<', np.max(np.fabs(si)), '|phi_i|<', np.max(np.fabs(phi_i))

    print '|phi_i|<', np.max(np.fabs(phi_i))

    if False:
        #pl.hist(si.flatten(), bins=20)
        pl.hist(phi_i.flatten(), bins=20, histtype='step')
        pl.show()


    return phi, phi_i



def shifted_ZA(p, si):
    ''' ->> shift uniform grid after obtaining '''

    _do_linear_grid_ = False

    if _do_linear_grid_:
        lgrid=np.linspace(0, p.boxsize, p.nbin, endpoint=False)
        grid=mar.meshgrid(lgrid, lgrid, lgrid)[::-1,...]
        print 'lgrid shape:', lgrid.shape, 'meshgrid shape:', grid.shape, si.shape
    else:
        _grid_=np.random.uniform(0., p.boxsize, 3*p.nbin)

        grid=mar.meshgrid(_grid_[:p.nbin], _grid_[p.nbin:2*p.nbin], _grid_[2*p.nbin:])[::-1,...]
        #print 'lgrid shape:', gridx.shape, 'meshgrid shape:', grid.shape, si.shape


    shifted=np.zeros((3, p.nbin, p.nbin, p.nbin))

    if displace_interpolation==False:

        for i in range(3):
            shift_=grid[i]+si[i]
        
            # ->> periodic boundary <<- #
            shift_[np.where(shift_<0.)] = p.boxsize+shift_[np.where(shift_<0.)]
            shift_[np.where(shift_>p.boxsize)] = shift_[np.where(shift_>p.boxsize)]-p.boxsize
            shifted[i]=np.copy(shift_)
        
            print 'shifting axis-', i, 'is done.'
    else:

        print 'si shape', si.shape, 
        shifted=pmv.particle_move(p,np.swapaxes(grid.reshape(3,p.nbin**3),0,1),si)
        shifted[np.where(shifted<0.)] = p.boxsize+shifted[np.where(shifted<0.)]
        shifted[np.where(shifted>p.boxsize)]=shifted[np.where(shifted>p.boxsize)]-p.boxsize

    if True:
        shift_=shifted.reshape(p.nbin,p.nbin,p.nbin, 3)
        pl.plot(shift_[:,:,100,1], shift_[:,:,100, 2], 'k.')

        #pl.plot(si[1,:,:,100], si[2,:,:,100], 'k.')
        #pl.plot(grid[1,:,:,100], grid[2,:,:,100], 'r.')
	pl.show()

    return shifted




def displaced_ZA(p, si, mpart):
    # -> reconstruction <<- #

    displaced=np.zeros((3, p.nbin, p.nbin, p.nbin))

    if displace_interpolation==False:

        for i in range(3):
            dis_=mpart[i]+si[i]

            # ->> periodic boundary <<- #
	    dis_[np.where(dis_<0.)] = p.boxsize+dis_[np.where(dis_<0.)]
	    dis_[np.where(dis_>p.boxsize)] = dis_[np.where(dis_>p.boxsize)]-p.boxsize
	    displaced[i] = np.copy(dis_)

    	    print 'displacing axis-', i, 'is done.'


    else:

        print 'si shape', si.shape, 'particle pos shape:', mpart.shape
        displaced=pmv.particle_move(p, mpart, si)

        displaced[np.where(displaced<0.)] = p.boxsize+displaced[np.where(displaced<0.)]
        displaced[np.where(displaced>p.boxsize)] = displaced[np.where(displaced>p.boxsize)]-p.boxsize

    return displaced





''' ->> Density Reconstruction Wrapper <<- '''
_rect_type_list_=['ZA_displaced', 'ZA_displaced_shifted']

def lag_rec_ZA(p, mpart, dmap, smooth_R=None, smooth_type=None, rect_type='ZA_displaced'): 
    # -> reconstruction of initial state, return density map <<- #

    # ->> get displacement field <<- #
    phi, si=get_ZA_displacement(p, dmap, smooth_R=smooth_R, smooth_type=smooth_type)
    print '->> phi, si shape:', phi.shape, si.shape

    npt=p.nbin**3.

    if rect_type=='ZA_displaced_shifted':
        #->> get particles <<- #
        #pos_disp=np.copy(np.swapaxes(displaced_ZA(p,si,mpart).reshape(3, p.nbin**3),0,1))
	#pos_shift=np.copy(np.swapaxes(shifted_ZA(p, si).reshape(3, p.nbin**3),0,1))

        pos_disp=displaced_ZA(p,si,mpart)
	pos_shift=shifted_ZA(p, si)

	print 'particle shape:', mpart.shape, pos_disp.shape, pos_shift.shape

        #->> converting density map <<- #
	d_disp = mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos_disp, pmass=p.particle_mass)
	d_shift= mcic.cic(p.cp, npt, p.nbin, p.boxsize, pos_shift, pmass=p.particle_mass)


        _test_draw_ = True
	if _test_draw_:
            shift_=pos_shift.reshape(p.nbin,p.nbin,p.nbin, 3)

            nplt, ncol = 4, 2
            fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=7.,gap_size=1,return_figure=True)


            ax[0].imshow(d_shift[:,:,100])
            ax[1].plot(shift_[:,:,100,1], shift_[:,:,100, 2], 'k.')

	    ax[2].imshow(1.+d_disp[:,:,100], norm=colors.LogNorm())
	    ax[3].imshow(1.+dmap[:,:,100], norm=colors.LogNorm())


	    pl.show()



	d_rec=d_disp-d_shift

        return d_rec, d_disp, d_shift, pos_disp, pos_shift


    elif rect_type=='ZA_displaced':

        return 


    else:
        raise Exception













        """
        # ->> setup interpolation of si <<- #
        x_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
        y_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
        z_= np.linspace(0., p.boxsize, p.nbin, endpoint=False)
    
        x, y, z=mar.meshgrid(x_, y_, z_)
        data =[np.array([x, y, z, si[i]]) for i in range(3)]
    
        print 'interpolation data shape:', x_.shape, x.shape, data[0].shape
    
        sitp=[]
        for i in range(3):
            sitp.append(itp.griddata_interpolater(data[i], 3, method='linear', fill_value=0))
    
    
        # ->> move particles <<- #
        for i in range(3):
            dis_=mpart[i]+sitp[i]( (mpart[0], mpart[1], mpart[2]) )

            # ->> periodic boundary <<- #
	    dis_[np.where(dis_<0.)] = p.boxsize+dis_[np.where(dis_<0.)]
	    dis_[np.where(dis_>p.boxsize)] = dis_[np.where(dis_>p.boxsize)]-p.boxsize
	    displaced[i] = np.copy(dis_)

    	    print 'axis-', i, 'is done.'
        """
