cdef extern from "../c/gadget_read_snapshot.h":

    int load_snapshot(char *fname, int files)

    int reordering(void)

