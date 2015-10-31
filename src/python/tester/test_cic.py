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








if __name__=='__main__':

    root='/mnt/scratch-lustre/xwang/data/velinv/cmpc/'
    fname_part='/gadget/snap100Mpc128_z0'
    fname_field='/field/snap100Mpc128_z0.den'

    #->> load data <<- #
    fn_part=root+fname_part
    s=pn.load(fn_part)

    nbin=128
    
    # ->> convert from kpc/h to Mpc/h <<- #
    pos = (s['pos']/1.e3).reshape(nbin,nbin,nbin,3)
    print 'pos.shape', pos.shape, 'pos boundary:', pos[...,0].min(), pos[...,0].max()


    ''' ->> run CIC density estimation <<- '''
    d=mcic.cic(pos, nbin)
    print d.shape, d.min(), d.max()


    if True:
        fig=pl.figure(figsize=(20, 20))
        ax=fig.add_subplot(111)

	data=d[...,100]+1e-3
        ax.imshow(np.flipud(data), norm=colors.LogNorm(vmin=data.min(),vmax=data.max()) )

        fig.savefig('dtest.png')


