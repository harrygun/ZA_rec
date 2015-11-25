#ifndef _H_POISSON_
#define _H_POISSON_


void poisson_solver(double *d, double *phi, int ngrid, int smooth_type, double smooth_R);
void poisson_solver_float(float *d, float *phi, int ngrid, int smooth_type, double smooth_R);

#endif
