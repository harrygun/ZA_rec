  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"
  #include "poisson.h"

  #define pi M_PI

  #ifdef _OMP_
  #include <omp.h>
  #endif


void smooth_field(float *d, double boxsize, int ngrid, int smooth_type, 
                                                       double smooth_R)  {
  /* ->> smooth field with FFT <<- */
  long long dksize, dsize, l, m, n, i, j;
  float kx, ky, kz, ki[3], sin2x, sin2y, sin2z, W, kmin;
  float fac=1.0/(float)(ngrid*ngrid*ngrid);
  kmin=2.*pi/boxsize;

  printf("\n->> Smoothing field with FFT.\n");

  // ->> initialize OpenMP <<- //
  #ifdef _OMP_
  if(fftwf_init_threads()==0) abort();
  #endif


  // ->> initialization <<- //
  dsize=ngrid*ngrid*ngrid*sizeof(float);
  dksize=ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex);
  fftwf_complex *dk;

  dk=(fftwf_complex *)fftwf_malloc(dksize);

  // ->> multi-threads initialization <<- //
  #ifdef _OMP_
  fftwf_plan_with_nthreads(omp_get_max_threads());
  #endif

  fftwf_plan pforward, pbackward;
  
  pforward=fftwf_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
  fftwf_execute(pforward);  

  
  /* smooth with FFT */
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

        // ->> window function <<- //
        if(smooth_type==_TOPHAT_SMOOTH_) { W=1.; }
        else if(smooth_type==_GAUSSIAN_SMOOTH_) {
          W=exp(-(sin2x+sin2y+sin2z)*smooth_R*smooth_R/2.); 
	  }
        else { W=1.; }

        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens*W;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens*W;
	  
        // ->> do gradient <<- //
        if(do_grad==TRUE) {
	  for (i=0; i<3; i++) {
            ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), i, l, m, n)[0]=
             ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*(-ki[i]);
            ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), i, l, m, n)[1]=
             ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*ki[i];
            }
	  }

        // ->>  do Hessian matrix <<- //
        if(do_hess==TRUE) {
	  for (i=0; i<3; i++)
	    for (j=0; j<3; j++) 
	      for (cc=0; cc<2; cc++) { // 2 complex components //
                ArrayAccess5D_n5(dkij, 3, 3, ngrid, ngrid, (ngrid/2+1), i, j, l, m, n)[cc]=
	              ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[cc]*ki[i]*ki[j]*(-1.);
	        }
	  }

        }
  
  /* find the inverse FFT of phi */

  if(do_phi==TRUE){
    printf("do phi itself.\n");
    pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, phi, FFTW_ESTIMATE);
    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);
    }

  int rank, howmany, *ndim, idist, odist, istride, ostride, *inembed, *onembed;
  ndim=(int *)malloc(3*sizeof(int));
  ndim[0]=ngrid; ndim[1]=ngrid; ndim[2]=ngrid;

  if (do_hess==TRUE) {
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


  //->> 
  if (do_grad==TRUE) {
    printf("do Gradient vector of phi.\n");

    rank=3;
    howmany=3;
    idist=ngrid*ngrid*(ngrid/2+1);
    odist=ngrid*ngrid*ngrid;
    istride=1; ostride=1;
    inembed=NULL; onembed=NULL;

    pbackward_grad=fftwf_plan_many_dft_c2r(rank, ndim, howmany, dki, inembed, istride, idist, phi_i, onembed, ostride, odist, FFTW_ESTIMATE);
    
    fftwf_execute(pbackward_grad);
    fftwf_destroy_plan(pbackward_grad);
    }
  
  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++) {

        if (do_phi==TRUE)
          ArrayAccess3D(phi, ngrid, l, m, n)*=fac;
  
        if (do_hess==TRUE) {
	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
              ArrayAccess5D_n5(phi_ij, 3, 3, ngrid, ngrid, ngrid, i, j, l, m, n)*=fac;
	  }

	if (do_grad==TRUE) {
	  for(i=0; i<3; i++)
            ArrayAccess4D_n4(phi_i, 3, ngrid, ngrid, ngrid, i, l, m, n)*=fac;
	  }

	}


  fftwf_free(dk);
  free(ndim);
  if(do_grad==TRUE) {fftwf_free(dki);}
  if(do_hess==TRUE) {fftwf_free(dkij);}

  fftwf_destroy_plan(pforward);
  fftwf_cleanup();
 
  printf("->> FFT is Done <<- \n\n");

  return;
  }









