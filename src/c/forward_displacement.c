  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "parvar.h"
  #include "io.h"
  #include "cic.h"
  #include "misc.h"
  #include "poisson.h"
  #include "fourier.h"
  #include "backward_displacement.h"


  #ifdef _OMP_
  #include <omp.h>
  #endif






void forward_displacement_one_iteration(SimInfo *s, float *d, float *disp_input, 
                                         float *disp_output){
  // ->> forward modeling of displacement field for single iteration <<- //






  return;
  }









/* ->>  Forward modeling of displacement field <<- */
void displacement_forward(SimInfo *s, float *d, float *disp, int nit_max) {
  /* ->> obtain ZA+2LPT displacement field from density field <<- */ 

  int fft_return_type;
  long long ip, i, j, ngrid_tot;
  float *d2lpt, *phi1, *phi1_i, *phi1_ij, *phi, *phi_i, *phi_ij;
  //char *other_req="";

  ngrid_tot=s->ngrid_xyz[0]*s->ngrid_xyz[1]*s->ngrid_xyz[2];

  // ->> solve FFTW  to get 1st order phi and phi_i <<- //
  fft_return_type=_RETURN_HESSIAN_;
  phi1_ij=(float *)fftwf_malloc(sizeof(float)*ngrid_tot*9);

  poisson_solver_float(d, phi1, phi1_i, phi1_ij, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R, fft_return_type);

  // ->> phi^(2) <<- //
  phi=(float *)fftwf_malloc(sizeof(float)*ngrid_tot);
  // ->> smooth density field first <<- //
  smooth_field(d, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R);

  #ifdef _OMP_
  #pragma omp parallel for private(ip, i, j)
  #endif
  for(ip=0; ip<ngrid_tot; ip++) {

    for(i=0; i<3; i++)
      for(j=0; j>i; j++) {
        d[ip] += 3./7.*(ArrayAccess3D_n3(phi1_ij, 3, 3, ngrid_tot, i, i, ip)
	           * ArrayAccess3D_n3(phi1_ij, 3, 3, ngrid_tot, j, j, ip) 
                   -pow(ArrayAccess3D_n3(phi1_ij, 3, 3, ngrid_tot, i, j, ip), 2));
        }
    }

  //->> solve Poisson equation again <<- //
  fft_return_type=_RETURN_GRADIENT_;

  poisson_solver_float(d, phi, disp, phi_ij, s->boxsize, s->ngrid, s->smooth_type_flag, 0., fft_return_type);


  return;
  }





