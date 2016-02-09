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
  #include "readfile.h"


  #include "parvar.h"
  #include "io.h"
  #include "cic.h"
  #include "misc.h"
  #include "poisson.h"
  #include "fourier.h"
  #include "backward_displacement.h"

  #include "stat_model.h"


  #define pi M_PI

  #ifdef _OMP_
  #include <omp.h>
  #endif





void potential_curlfree_vec(float *disp, double boxsize, int ngrid){
  //->> return both the potential and the curl-free part of give vector field <<- //
  int rank, howmany, *ndim, idist, odist, istride, ostride, *inembed, *onembed;
  long long dksize, dsize, l, m, n, i, j;
  float kx, ky, kz, ki[3], sin2x, sin2y, sin2z, W, kmin;
  float fac=1.0/(float)(ngrid*ngrid*ngrid);
  kmin=2.*pi/boxsize;

  printf("\n->> Smoothing field with FFT.\n");

  // ->> initialize OpenMP <<- //
  #ifdef _OMP_
  if(fftwf_init_threads()==0) abort();
  #endif

  // ->>   <<- //
  // ->> initialization <<- //
  dsize=ngrid*ngrid*ngrid*sizeof(float);
  dksize=ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex);
  fftwf_complex *dki;

  dki=(fftwf_complex *)fftwf_malloc(3*dksize);

  // ->> multi-threads initialization <<- //
  #ifdef _OMP_
  fftwf_plan_with_nthreads(omp_get_max_threads());
  #endif

  fftwf_plan pforward, pbackward;
  
  rank=3;
  howmany=3;
  ndim=(int *)malloc(3*sizeof(int));
  ndim[0]=ngrid; ndim[1]=ngrid; ndim[2]=ngrid;

  idist=ngrid*ngrid*ngrid;
  odist=ngrid*ngrid*(ngrid/2+1);

  istride=1; ostride=1;
  inembed=NULL; onembed=NULL;


  pforward=fftwf_plan_many_dft_r2c(rank, ndim, howmany, dki, inembed, istride, idist, phi_ij, onembed, ostride, odist, FFTW_ESTIMATE);

  //pforward=fftwf_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk, FFTW_ESTIMATE);
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
	

        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=W;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=W;
        }
  
    /* inverse FFT */
    pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, d, FFTW_ESTIMATE);
    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);

  
  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++) 
        ArrayAccess3D(d, ngrid, l, m, n)*=fac;

  fftwf_free(dk);

  fftwf_destroy_plan(pforward);
  fftwf_cleanup();
 
  printf("->> FFT is Done <<- \n\n");
  return;
  }










void smooth_field(float *d, double boxsize, int ngrid, int smooth_type, 
                                          double smooth_R, Interpar *sw)  {
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

        // ->> window function <<- //
        if(smooth_type==_TOPHAT_SMOOTH_) { W=1.; }
        else if(smooth_type==_GAUSSIAN_SMOOTH_) {
          W=exp(-(sin2x+sin2y+sin2z)*smooth_R*smooth_R/2.); 
	  }
	else if(smooth_type==_INVERSE_GAUSSIAN_SMOOTH_){
          W=1./(1.-exp(-(sin2x+sin2y+sin2z)*smooth_R*smooth_R/2.)); 
	  }
	else if(smooth_type==_ANISOTROPIC_INTERPOLATOR_SMOOTH_) {
	  W=tk_interp(sw, sqrt(sin2x+sin2y+sin2z));
          //printf("W[%lg, %lg, %lg]=%lg\n", ki[0], ki[1], ki[2], W);
	  }
        else { 
          printf("smooth_field window function error.\n");  fflush(stdout);
	  W=1.; }

        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=W;
        ArrayAccess3D_n3(dk, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=W;
        }
  
    /* inverse FFT */
    pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk, d, FFTW_ESTIMATE);
    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);

  
  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++) 
        ArrayAccess3D(d, ngrid, l, m, n)*=fac;

  fftwf_free(dk);

  fftwf_destroy_plan(pforward);
  fftwf_cleanup();
 
  printf("->> FFT is Done <<- \n\n");
  return;
  }



