  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "parvar.h"

  #define pi M_PI



void poisson_solver(double *d, double *phi, int ngrid, int smooth_type, double smooth_R)  {
  /* ->> Poisson Solver with FFT <<- */

  int l, m, n;
  double kx, ky, kz, sin2x, sin2y, sin2z, greens, W;
  double fac=1.0/(double)(ngrid*ngrid*ngrid);

  fftw_complex *dk=(fftw_complex *)fftw_malloc(ngrid*ngrid*(ngrid/2+1)*sizeof(fftw_complex));

  fftw_plan pforward, pbackward;
  
  pforward=fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
  fftw_execute(pforward);  
  
  /* work out the green's function and the FFT of phi*/
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid/2+1; n++){
        kx = 2.*pi*l/(double)ngrid;
        ky = 2.*pi*m/(double)ngrid;
        kz = 2.*pi*n/(double)ngrid;
        
        sin2x = sin(kx/2.)*sin(kx/2.);
        sin2y = sin(ky/2.)*sin(ky/2.);
        sin2z = sin(kz/2.)*sin(kz/2.);
        
        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = -1./(sin2x+sin2y+sin2z);
      	    
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens;
        }
  
  /* find the inverse FFT of phi */
  pbackward = fftw_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, phi, FFTW_ESTIMATE);
  fftw_execute(pbackward);
  
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++)
        ArrayAccess3D(phi, ngrid, l, m, n)*=fac;
  
  fftw_destroy_plan(pforward);
  fftw_destroy_plan(pbackward);

  return;
  }




void poisson_solver_float(float *d, float *phi, int ngrid, int smooth_type, double smooth_R)  {
  /* ->> Poisson Solver with FFT <<- */

  int l, m, n;
  double kx, ky, kz, sin2x, sin2y, sin2z, greens, W;
  double fac=1.0/(double)(ngrid*ngrid*ngrid);

  fftwf_complex *dk=(fftwf_complex *)fftwf_malloc(ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex));

  fftwf_plan pforward, pbackward;
  
  pforward=fftwf_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
  fftwf_execute(pforward);  
  
  /* work out the green's function and the FFT of phi*/
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid/2+1; n++){
        kx = 2.*pi*l/(double)ngrid;
        ky = 2.*pi*m/(double)ngrid;
        kz = 2.*pi*n/(double)ngrid;
        
        sin2x = sin(kx/2.)*sin(kx/2.);
        sin2y = sin(ky/2.)*sin(ky/2.);
        sin2z = sin(kz/2.)*sin(kz/2.);
        
        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = -1./(sin2x+sin2y+sin2z);
      	    
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens;
        }
  
  /* find the inverse FFT of phi */
  pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, phi, FFTW_ESTIMATE);
  fftwf_execute(pbackward);
  
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++)
        ArrayAccess3D(phi, ngrid, l, m, n)*=fac;
  
  fftwf_destroy_plan(pforward);
  fftwf_destroy_plan(pbackward);

  return;
  }
