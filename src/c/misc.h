#ifndef _H_MISC_
#define _H_MISC_


void mat_inv_3d_float(float mat[3][3], float imat[3][3]);

void mat_multiply_3d_float(float mat[3][3], float x[3], float y[3]);


void cubic_trim(float *mi, float *mo, size_t ngrid, size_t ntrim);

#endif
