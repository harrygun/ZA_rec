import numpy as np
import matplotlib.pyplot as pl
import pynbody as pn

import genscript.read as rd
import io as mio






'''-------------------------------------------------------------------------- '''

def import_MIP_data(p, droot_part, droot_field, fn, import_data_type='all'):

    print ' >> Importing data ...'

    if import_data_type=='all':
        datype=['particle', 'field']
    else:
        datype=[import_data_type]

    if 'particle' in datype:
        npart=int(p.nbin**3)
        _x, _y, _z = mio.read_cgal(droot_part, fn, npart, import_type='position')
        print _x.shape

        x=_x.reshape(p.nbin, p.nbin, p.nbin)*p.boxsize
        y=_y.reshape(p.nbin, p.nbin, p.nbin)*p.boxsize
        z=_z.reshape(p.nbin, p.nbin, p.nbin)*p.boxsize
        pos=np.array([x, y, z])
    
    # ->> import field <<- #
    if 'field' in datype:
        fn_d = droot_field+fn+'.'+str(p.nbin)+'.fvol'
        den=rd.rgrid(fn_d, ngrid=p.nbin, dtype='float', comp=1)
        print 'imported density shape:', den.shape

    if import_data_type=='all':
        return pos, den
    elif import_data_type=='particle':
        return pos
    elif import_data_type=='field':
        return den


def import_gadget_DTFE(p, fn_part, fn_field, import_data_type='all'):

    s=pn.load(fn_part)
    # ->> convert from kpc/h to Mpc/h <<- #
    pos = np.reshape(np.swapaxes(s['pos']/1.e3, 0, 1), (3, p.nbin, p.nbin, p.nbin))
    print 'pos.shape', pos.shape, 'pos boundary:', pos[0].min(), pos[1].max()

    # ->> field
    den=rd.rgrid(fn_field, ngrid=p.nbin, dtype='float', comp=1)
    print 'den shape:', den.shape

    if import_data_type=='all':
        return pos, den
    elif import_data_type=='particle':
        return pos
    elif import_data_type=='field':
        return den



def import_cita_simulation(p, fn_part, fn_field, import_data_type='all'):

    if import_data_type in ['all', 'field']:
        _import_field_=True
    else:
        _import_field_=False

    pos_, v=mio.read_cita_simulation(fn_part, p.nbin)

    xmin, xmax =np.zeros(3), np.zeros(3)
    pos=np.zeros(pos_.shape)

    '''
    for i in range(3):
        xmax[i], xmin[i] = np.max(pos_[...,i]), np.min(pos_[...,i])
	pos[...,i]=(pos_[...,i]-xmin[i])*p.boxsize/(xmax[i]-xmin[i]) 
    '''
    pos=pos_

    #pos=np.rollaxis(pos, -1).reshape(3, p.nbin, p.nbin, p.nbin)

    #print 'shape:', den.shape, pos.shape
    print 'pos max/min:', xmax, xmin


    if _import_field_:
        den=np.load(fn_field)['d']
	return  pos, den
    else:
	return  pos, None


