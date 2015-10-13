import numpy as np
cimport numpy as np
cimport gadget_read as gr




def read_snapshot(fname, nfiles):
    ''' `files`: a snapshot may be divided into `nfiles` files.  '''

    gr.load_snapshot(<char *>fname, <int> files)

    return 
