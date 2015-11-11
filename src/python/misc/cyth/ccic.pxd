cdef extern from "../c/ccic.h":

    double density(float *Pos, float *delta, long long NDM, long long Ngrid, double mass) 

