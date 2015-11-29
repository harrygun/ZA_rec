  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"
  #include "poisson.h"

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

  int l, m, n, i, j, cc, do_grad, do_hess, dksize, dsize;
  float kx, ky, kz, ki[3], sin2x, sin2y, sin2z, greens, W, kmin;
  float fac=1.0/(float)(ngrid*ngrid*ngrid);
  kmin=2.*pi/boxsize;

  // ->> return type <<- //
  if((return_type==_RETURN_GRADIENT_)||(return_type==_RETURN_GRADIENT_HESSIAN_))
    do_grad=TRUE;
  else 
    do_grad=FALSE;
  if((return_type==_RETURN_HESSIAN_)||(return_type==_RETURN_GRADIENT_HESSIAN_))
    do_hess=TRUE;
  else 
    do_hess=FALSE;

  // ->> initialization <<- //
  dsize=ngrid*ngrid*ngrid*sizeof(float);
  dksize=ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex);
  fftwf_complex *dk, *dki, *dkij;

  dk=(fftwf_complex *)fftwf_malloc(dksize);

  if(do_grad==TRUE)
    dki=(fftwf_complex *)fftwf_malloc(3*dksize);
  if(do_hess==TRUE)
    dkij=(fftwf_complex *)fftwf_malloc(3*3*dksize);

  fftwf_plan pforward, pbackward, pbackward_hess;
  
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
	ki[0]=2.*sin(kx/2.);
	ki[1]=2.*sin(ky/2.);
	ki[2]=2.*sin(kz/2.);
	
        //sin2x = kx*kx; sin2y = ky*ky; sin2z = kz*kz;
	//ki[0]=kx; ki[1]=ky; ki[2]=kz;

        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = -1./(sin2x+sin2y+sin2z);
      	    
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens;

        // ->> do gradient <<- //
        if(do_grad) {
	  for (i=0; i<3; i++)
            for (cc=0; cc<2; cc++) { 
                ArrayAccess4D_n4(dkij, 3, ngrid, ngrid, (ngrid/2+1), i, j, l, m, n)[cc]=
	              ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[cc]*ki[i]*ki[j]*(-1.);

	      }
	  }

        // ->>  do Hessian matrix <<- //
        if(do_hess) {
	  for (i=0; i<3; i++)
	    for (j=0; j<3; j++) 
	      for (cc=0; cc<2; cc++) { // 2 complex components //
                ArrayAccess5D_n5(dkij, 3, 3, ngrid, ngrid, (ngrid/2+1), i, j, l, m, n)[cc]=
	              ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[cc]*ki[i]*ki[j]*(-1.);
	        }
	  }

        }
  
  /* find the inverse FFT of phi */
  pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, phi, FFTW_ESTIMATE);
  fftwf_execute(pbackward);

  int rank, howmany, *ndim, idist, odist, istride, ostride, *inembed, *onembed;
  ndim=(int *)malloc(3*sizeof(int));
  ndim[0]=ngrid; ndim[1]=ngrid; ndim[2]=ngrid;

  if (do_hess) {
    printf("do Hessian matrix of phi.\n");

    rank=3;
    howmany=9;
    idist=ngrid*ngrid*(ngrid/2+1);
    odist=ngrid*ngrid*ngrid;
    istride=1; ostride=1;
    inembed=NULL; onembed=NULL;

    pbackward_hess=fftwf_plan_many_dft_c2r(rank, ndim, howmany, dkij, inembed, istride, idist, phi_ij, onembed, ostride, odist, FFTW_ESTIMATE);
    
    fftwf_execute(pbackward_hess);
    fftwf_destroy_plan(pbackward_hess);
    }

  
  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++) {
        ArrayAccess3D(phi, ngrid, l, m, n)*=fac;
  
        if (do_hess) {
	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
              ArrayAccess5D_n5(phi_ij, 3, 3, ngrid, ngrid, ngrid, i, j, l, m, n)*=fac;
	  }

	if (do_grad) {
	  }

	}


  fftwf_free(dk);
  if(do_grad) fftwf_free(dki); 
  if(do_hess) fftwf_free(dkij);

  fftwf_destroy_plan(pforward);
  fftwf_destroy_plan(pbackward);

  //fftwf_cleanup();

  return;
  }




