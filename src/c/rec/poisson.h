#ifndef _H_POISSON_
#define _H_POISSON_


  #define _RETURN_PHI_ONLY_  2121
  #define _RETURN_GRADIENT_  2121
  #define _RETURN_HESSIAN_   2122
  #define _RETURN_GRADIENT_HESSIAN_  2123

//void poisson_solver(double *d, double *phi, double boxsize, int ngrid, int smooth_type, double smooth_R); 
//void poisson_solver_float(float *d, float *phi, double boxsize, int ngrid, int smooth_type, double smooth_R);



void poisson_solver_float(float *d, float *phi, float *phi_i, float *phi_ij, double boxsize, int ngrid,
                           int smooth_type, double smooth_R, int return_type);
#endif
