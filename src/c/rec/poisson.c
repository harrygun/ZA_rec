  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "parvar.h"

  #define pi M_PI



void poisson_solver(double *d, double *phi, double boxsize, int ngrid, int smooth_type, double smooth_R)  {
  /* ->> Poisson Solver with FFT <<- */

  int l, m, n;
  double kx, ky, kz, sin2x, sin2y, sin2z, greens, W, kmin;
  double fac=1.0/(double)(ngrid*ngrid*ngrid);

  kmin=2.*pi/boxsize;

  fftw_complex *dk=(fftw_complex *)fftw_malloc(ngrid*ngrid*(ngrid/2+1)*sizeof(fftw_complex));

  fftw_plan pforward, pbackward;
  
  pforward=fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
  fftw_execute(pforward);  
  
  /* work out the green's function and the FFT of phi*/
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid/2+1; n++){

        if(l<ngrid/2) kx=l*kmin;
	else kx=(l-ngrid)*kmin;

        if(m<ngrid/2) ky=m*kmin;
	else ky=(m-ngrid)*kmin;

        if(n<ngrid/2) kz=n*kmin;
	else kz=(n-ngrid)*kmin;

        
        sin2x = 4.*sin(kx/2.)*sin(kx/2.);
        sin2y = 4.*sin(ky/2.)*sin(ky/2.);
        sin2z = 4.*sin(kz/2.)*sin(kz/2.);
        
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
  
  fftw_free(dk);
  fftw_destroy_plan(pforward);
  fftw_destroy_plan(pbackward);

  return;
  }




void poisson_solver_float(float *d, float *phi, float *phi_i, float *phi_ij, double boxsize, int ngrid,
                           int smooth_type, double smooth_R, int return_type)  {
  /* ->> Poisson Solver with FFT <<- */

  int l, m, n, i, j, do_grad, do_hess;
  double kx, ky, kz, _kx, _ky, _kz, sin2x, sin2y, sin2z, greens, W, kmin;
  double fac=1.0/(double)(ngrid*ngrid*ngrid);
  kmin=2.*pi/boxsize;

  // ->> return type <<- //
  if((return_type==_RETURN_GRADIENT_)&&(return_type==_RETURN_GRADIENT_HESSIAN_))
    do_grad=TRUE;
  else 
    do_grad=FALSE
  if((return_type==_RETURN_HESSIAN_)&&(return_type==_RETURN_GRADIENT_HESSIAN_))
    do_hess=TRUE;
  else 
    do_hess=FALSE

  // ->> initialization <<- //
  fftwf_complex *dk=(fftwf_complex *)fftwf_malloc(ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex));
  if(do_grad==TRUE)
    fftwf_complex *dki=(fftwf_complex *)fftwf_malloc(ngrid*ngrid*(ngrid/2+1)*3*sizeof(fftwf_complex));
  if(do_hess==TRUE)
    fftwf_complex *dkij=(fftwf_complex *)fftwf_malloc(ngrid*ngrid*(ngrid/2+1)*3*3*sizeof(fftwf_complex));

  fftwf_plan pforward, pbackward;
  
  pforward=fftwf_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
  fftwf_execute(pforward);  
  
  /* work out the green's function and the FFT of phi*/
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid/2+1; n++){

        if(l<ngrid/2) kx=l*kmin;
	else kx=(l-ngrid)*kmin;

        if(m<ngrid/2) ky=m*kmin;
	else ky=(m-ngrid)*kmin;

        if(n<ngrid/2) kz=n*kmin;
	else kz=(n-ngrid)*kmin;

        sin2x = 4.*sin(kx/2.)*sin(kx/2.);
        sin2y = 4.*sin(ky/2.)*sin(ky/2.);
        sin2z = 4.*sin(kz/2.)*sin(kz/2.);

        //sin2x = kx*kx;
        //sin2y = ky*ky;
        //sin2z = kz*kz;
        
        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = -1./(sin2x+sin2y+sin2z);
      	    
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens;

        if(do_grad) {
	  for (i=0; i<3; i++)
	  }

        if(do_hess) {
	  }

        }
  
  /* find the inverse FFT of phi */
  pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, phi, FFTW_ESTIMATE);
  fftwf_execute(pbackward);
  
  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++)
        ArrayAccess3D(phi, ngrid, l, m, n)*=fac;


  fftwf_free(dk);
  if(do_grad) fftwf_free(dki); 
  if(do_hess) fftwf_free(dkij);

  fftwf_destroy_plan(pforward);
  fftwf_destroy_plan(pbackward);

  return;
  }




