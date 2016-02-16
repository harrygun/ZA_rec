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





void potential_curlfree_vec(float *disp, float *div, float *phi, float *disp_phi, double boxsize, int ngrid) {
  /*->> return both the potential and the curl-free part of give vector field <<- */
  // ->> disp: (input)  original vector field, presumably displacement field 
  // ->> div:  (output) the divergence of disp 
  // ->> phi:  (output) the potential field of disp so that \nabla_i phi= disp_i
  // ->> disp_phi: (output)  the curl-free part of disp 

  int do_phi, do_div, do_disp_phi;

  do_phi=TRUE; do_div=TRUE, do_disp_phi=TRUE;
  if(div==NULL) { do_div=FALSE; }
  if(phi==NULL) { do_phi=FALSE; }
  if(disp_phi==NULL) { do_disp_phi=FALSE; }
  
  if( (do_div==FALSE)&&(do_phi==FALSE)&&(do_disp_phi==FALSE) )
    abort();

  int rank, howmany, *ndim, idist, odist, istride, ostride, *inembed, *onembed;
  long long dksize, dsize, l, m, n, i, j;
  float kx, ky, kz, ki[3], sin2x, sin2y, sin2z, W, kmin, greens;
  float fac=1.0/(float)(ngrid*ngrid*ngrid);
  kmin=2.*pi/boxsize;

  printf("\n->> Obtain potential field with FFT.\n");

  // ->> initialize OpenMP <<- //
  #ifdef _OMP_
  if(fftwf_init_threads()==0) abort();
  #endif

  // ->>   <<- //
  // ->> initialization <<- //
  dsize=ngrid*ngrid*ngrid*sizeof(float);
  dksize=ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex);
  fftwf_complex *dki, *dkdiv, *dkphi;

  dki=(fftwf_complex *)fftwf_malloc(3*dksize);
  dkdiv=(fftwf_complex *)fftwf_malloc(dksize);
  dkphi=(fftwf_complex *)fftwf_malloc(dksize);

  // ->> multi-threads initialization <<- //
  #ifdef _OMP_
  fftwf_plan_with_nthreads(omp_get_max_threads());
  #endif

  //->> forward FFT <<- //
  fftwf_plan pforward, pbackward;
  
  rank=3; howmany=3;

  ndim=(int *)malloc(3*sizeof(int));
  ndim[0]=ngrid; ndim[1]=ngrid; ndim[2]=ngrid;

  idist=ngrid*ngrid*ngrid;
  odist=ngrid*ngrid*(ngrid/2+1);

  istride=1; ostride=1;
  inembed=NULL; onembed=NULL;

  pforward=fftwf_plan_many_dft_r2c(rank, ndim, howmany, disp, inembed, istride, idist, dki, onembed, ostride, odist, FFTW_ESTIMATE);

  fftwf_execute(pforward);  
  //->> end of forward FFT <<-//

  
  /* ->> taking the divergence, potential, and curl-free vector <<- */
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

        //->> divergence <<- //
        ArrayAccess3D_n3(dkdiv, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]=
               -ki[0]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 0, l, m, n)[1]+
               -ki[1]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 1, l, m, n)[1]+
               -ki[2]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 2, l, m, n)[1];

        ArrayAccess3D_n3(dkdiv, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]=
                ki[0]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 0, l, m, n)[0]+
                ki[1]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 1, l, m, n)[0]+
                ki[2]*ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), 2, l, m, n)[0];

        // ->> potential <<- //
        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = -1./(sin2x+sin2y+sin2z);

        ArrayAccess3D_n3(dkphi, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]=
                ArrayAccess3D_n3(dkdiv, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*greens;
        ArrayAccess3D_n3(dkphi, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]=
                ArrayAccess3D_n3(dkdiv, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*greens;


	//->> get the curl-free vector <<- //
	for(i=0; i<3; i++) {
          ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), i, l, m, n)[0]=
                 -ki[i]*ArrayAccess3D_n3(dkphi, ngrid, ngrid, (ngrid/2+1), l, m, n)[1];
          ArrayAccess4D_n4(dki, 3, ngrid, ngrid, (ngrid/2+1), i, l, m, n)[1]=
                  ki[i]*ArrayAccess3D_n3(dkphi, ngrid, ngrid, (ngrid/2+1), l, m, n)[0];
	  }

        }
  
  /* ->> inverse FFT <<- */

  // ->> get divergence first <<- //
  if(do_div==TRUE) {
    pbackward= fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dkdiv, div, FFTW_ESTIMATE);
    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);
    }

  // ->> then the potential phi <<- //
  if(do_phi==TRUE) {
    pbackward = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dkphi, phi, FFTW_ESTIMATE);
    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);
    }

  // ->> finally the curl-free vector field <<- //
  if(do_disp_phi==TRUE) {
    idist=ngrid*ngrid*(ngrid/2+1);
    odist=ngrid*ngrid*ngrid;

    pbackward=fftwf_plan_many_dft_c2r(rank, ndim, howmany, dki, inembed, istride, idist, disp_phi, onembed, ostride, odist, FFTW_ESTIMATE);

    fftwf_execute(pbackward);
    fftwf_destroy_plan(pbackward);
    }


  // ->> renormalize <<- //
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++) {

        if(do_div==TRUE) {ArrayAccess3D(div, ngrid, l, m, n)*=fac;}
        if(do_phi==TRUE) {ArrayAccess3D(phi, ngrid, l, m, n)*=fac;}

        if(do_disp_phi==TRUE) {
          for(i=0; i<3; i++)
            ArrayAccess4D_n4(disp_phi, 3, ngrid, ngrid, ngrid, i, l, m, n)*=fac;
	  }
	}



  // ->> free <<- //
  fftwf_free(dki); fftwf_free(dkdiv); fftwf_free(dkphi);
  free(ndim);

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
      for (n=0; n<ngrid; n++) {
        ArrayAccess3D(d, ngrid, l, m, n)*=fac;
	}

  fftwf_free(dk);

  fftwf_destroy_plan(pforward);
  fftwf_cleanup();
 
  printf("->> FFT is Done <<- \n\n");
  return;
  }



