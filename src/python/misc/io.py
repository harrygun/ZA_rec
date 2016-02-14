import os, sys
import array as arr
import numpy as np
import cyth.cgalio as cg
import genscript.mpiutil as mpi
from genscript.extendclass import *
import genscript.progcontrol as pc

#import pynbody as pn



def fnames(folder):
    files=os.listdir(folder)
    prefix=[ files[i].split('_') for i in range(len(files)) ]
    idx=np.where([ prefix[i][0]=='run' for i in range(len(prefix))])[0]
    f=np.sort(np.array(files)[idx])

    return f


def read_cgal(root, fn, npart, import_type='position'):

    if import_type=='all':
        _type_=['position', 'velocity']
    else:
        _type_=[import_type]

    if 'position' in _type_:
        fname=root+fn+'.CGAL'
        x, y, z= cg.read_position(fname, npart)

    if 'velocity' in _type_:
        sf=['.vx', '.vy', '.vz']
        vd=np.zeros((3, npart)).astype(np.float64)
        for i in range(3):
            fname=root+fn+sf[i]
            vd[i]=cg.read_velocity(fname, npart)

    if import_type=='all':
        return x, y, z, vd
    elif import_type=='position':
        return x, y, z
    elif import_type=='velocity':
        return vd




def read_cita_simulation(fn, npt):
    print 'importing data...'
    F=open(fn, 'rb')

    ''' ->> indicating the header <<- 
    header=arr.array('i')
    header.fromfile(F, 1)
    check=arr.array('f')
    check.fromfile(F, 11)
    print 'header:', header, header[0]
    '''

    numpy_type= np.float32

    head=F.read(48)
    #d = arr.array('f')
    #d.fromfile(F, npt**3*6)

    d=np.fromfile(F, dtype=np.dtype(numpy_type)).astype(numpy_type)

    #d=np.array(d).reshape(npt**3,6)
    #pos, vel=d[...,:3], d[...,3:6]

    F.close()

    print 'file closed.'

    #return pos, vel
    return np.array(d).reshape(npt**3,6)[...,:3], np.array(d).reshape(npt**3,6)[...,3:6]


def read_cita_simulation_pid(fn, ngrid):
    print 'importing PID data ..., ngrid=', ngrid

    F=open(fn, 'rb')
    numpy_type=np.int64

    #head=F.read(48)
    pid=np.fromfile(F, dtype=np.dtype(numpy_type), count=ngrid**3).astype(numpy_type)
    F.close()

    print 'PID file closed.', len(pid), ngrid**3

    return pid




"""
def gadget2cgal(root, fn):

    # ->> particle positions
    fname=root+fn

    s=pn.load(fname)
    pos = s['pos'].view(type=np.ndarray)
    v   = s['vel'].view(type=np.ndarray)

    # ->> write positions
    output=root+fn+'.CGAL'
    cg.write_position(output, pos[:,0], pos[:,1], pos[:,2])


    # ->> particle velocity
    sf=['.vx', '.vy', '.vz']
    for i in range(3):
	output=root+fn+sf[i]
	cg.write_velocity(output, v[:,i])


    return 
"""


def convert_gadget_all(p, root):
    f=fnames(root)
    print f

    if p.mpi:
        # >> MPI parallelization
	frange=mpi.mpirange(len(f))
    else:
        frange=range(len(f))

    for i in frange:
        print i, f[i]
        gadget2cgal(root, f[i])

    return





def cgal2delaunay(p, root, files, nbin):
    dtfe='/home/wangxin/workspace/code/velocity/workspace/data/Miguel/MPI-sim/cgal_dtfe_fast_vels '

    if p.mpi:
        # >> MPI parallelization
	frange=mpi.mpirange(len(files))
    else:
        frange=range(len(files))

    nbs=' '+str(nbin)+' '
    

    for i in frange:
        print i, files[i] 
        cmd = dtfe+root+files[i]+nbs+nbs+nbs+nbs+' 0 0 0 0.01'

	print cmd
        #os.system(cmd)


    return







if __name__=='__main__':

    init_dict=myDict({})
    p=pc.prog_init(**init_dict)


    root='/home/wangxin/workspace/code/velocity/workspace/data/Miguel/MPI-sim/'

    if False:
        # convert all gadget file to CGAL file
        convert_gadget_all(p, root)

    
    if False:
        # >> compute Delaunay tessellation from CGAL files 
        files=fnames(root)
        root='/home/wangxin/workspace/code/velocity/workspace/data/Miguel/MPI-sim/CGAL/'
        cgal2delaunay(p, root, files, 128)


    if True:
        # ->>
	#root='/home/wangxin/workspace/code/velocity/workspace/data/zeldovich/delaunay/256/gadget/'
	#fn='ics_z100'  #'ics_z2_256'
	root='/datascope/velinv/evolution/gadget/'
	fn='run_001'

        gadget2cgal(root, fn)




    if False:
        # ->> test 
	root='/home/wangxin/workspace/code/velocity/workspace/data/zeldovich/delaunay/256/gadget/'
	fn='ics_z100'  #'ics_z2_256'

        s=pn.load(root+fn)
	pos=s['pos']
	vel=s['vel']
        
	# ->> writing
        output=root+fn+'__test.CGAL'
        cg.write_position(output, pos[:,0], pos[:,1], pos[:,2])

        npart, boxsize=cg.read_header(output)
	print 'npart=', npart, 'boxsize=', boxsize

        # ->> reading
        x, y, z= cg.read_position(output, npart)


        vv=np.zeros((3, npart)).astype(np.float64)
        for i in range(3):
	    # ->> writing
            output=root+fn+'__test.v'+str(i+1)
            cg.write_velocity(output, vel[:,i])

            # ->> reading
	    vv[i]=cg.read_velocity(output, npart)


	print pos[:,1]
	print y

        print '=velocity='
	print vel[:,1]
	print vv[1,:]


    p.finalize()
