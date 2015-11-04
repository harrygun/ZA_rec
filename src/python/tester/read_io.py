import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
import matplotlib.colors as colors

import numpy as np
import array as arr
import pynbody as pn

import misc.io as mio
import misc.cyth.cgalio as cg
import misc.cic as mcic


def read_test(fn, npt):

    F=open(fn, 'rb')

    head=F.read(48)
    _d = arr.array('f')
    _d.fromfile(F, npt**3*6)

    #d=np.array(_d).reshape(npt,npt,npt,6)
    d=np.array(_d).reshape(npt**3,6)
    pos, vel=d[...,:3], d[...,3:6]

    F.close()

    return pos, vel


def write_cgal(fn, pos, vel):

    # ->> position <<- #
    output=fn+'.CGAL'
    cg.write_position(output, pos[:,0], pos[:,1], pos[:,2])

    # ->> velocity <<- #
    sf=['.vx', '.vy', '.vz']
    for i in range(3):
	output=fn+sf[i]
	cg.write_velocity(output, vel[:,i])

    return


def write_gadget(fn, pos):
    # ->> write to gadget file <<- #

    return



if __name__=='__main__':

    root='/mnt/scratch-lustre/xwang/data/baorec/cubep3m_dm_sml/node0/'
    fn_import='0.000xv0.dat'

    #fn='delta0'
    fn=root+fn_import
    npt=256

    #pos, vel =read_test(fn,npt)
    pos, vel = mio.read_cita_simulation(fn, npt)
    print 'pos shape', pos.shape


    #->>  <<- #
    if False:
        pos=pos.reshape(npt,npt,npt,3)

        fig=pl.figure(figsize=(20, 20))
        ax=fig.add_subplot(111)
        ax.plot(pos[:,:,100,1], pos[:,:,100,2], 'k.')
        #pl.show()
        fig.savefig('test.png')


    if False:
        # ->> convert into CGAL format <<- #
        fn=root+fn_import
        write_cgal(fn, pos, vel)

    if False:
        s=pn.load(fn_part)
        # ->> convert from kpc/h to Mpc/h <<- #
        pos = s['pos']/1.e3
        print 'pos.shape', pos.shape, 'pos boundary:', pos[0].min(), pos[1].max()


    if True:
        # ->> testing CIC <<- #
        npart=pos.shape[0]
	nbin=npt
        d=mcic.cic(npart, nbin, pos, pmass=1.e5)
	print 'd shape:', d.shape

        fn_out='0.000xv0.dat.den.npz'
	np.savez(fn_out, d=d)
         
        if False:
            fig=pl.figure(figsize=(20, 20))
            ax=fig.add_subplot(111)

	    data=d[...,100]+1e-3
            ax.imshow(np.flipud(data), norm=colors.LogNorm(vmin=data.min(),vmax=data.max()) )

            fig.savefig('cita_test.png')
