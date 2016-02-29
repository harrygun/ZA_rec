  #include <stdio.h>
  #include <time.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <gsl/gsl_histogram2d.h>
  #include <iniparser.h>

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
  #include "reconstruction_partmoving.h"
  #include "test_likelihood.h"

  #include "stat_model.h"
  #include "maxlike_phi_rec.h"
  #include "fourier.h"


#ifdef _MPI_
  #include "mpi.h"
#endif

#ifdef _OMP_
  #include <omp.h>
#endif





/* ->> max likelihood curve <<- */
int phi_mlik_init(Interpar *mlik, char *fname){

  FILE *fp;
  if(!(fp=fopen(fname, "r"))) {
    return FALSE; }

  // ->> importing data <<- //
  int i, j, line, extdidx=10;
  line=countFILEline(fp);

  double *phi, *mlik_phi, *smlik_phi, p, dp, dm;
  phi=(double *)malloc(line*sizeof(double));
  mlik_phi=(double *)malloc(line*sizeof(double));
  smlik_phi=(double *)malloc(line*sizeof(double));

  for(i=0; i<line; i++) 
    fscanf(fp, "%lg  %lg  %lg\n", &phi[i], &mlik_phi[i], &smlik_phi[i]);


  // ->> initialize interpolator <<- //
  myinterp_init(mlik, phi, smlik_phi, line);

  //->> linear extrapolation initialization <<- //
  mlik->min=phi[0];
  mlik->max=phi[line-1];

  dp=phi[extdidx]-phi[0]; dm=smlik_phi[extdidx]-smlik_phi[0];
  mlik->slop_min=dm/dp;

  dp=phi[line-extdidx]-phi[line-1]; dm=smlik_phi[line-extdidx]-smlik_phi[line-1];
  mlik->slop_max=dm/dp;

  //#define _DO_MLIK_TEST_
  #ifdef _DO_MLIK_TEST_
  fp=fopen("result/test_mlik.dat", "w");

  for(i=0; i<500; i++) {
    p=-300+i*600./500.; 
    printf("%lg  %lg\n", p, mlik_interp(mlik, p) );
    fprintf(fp, "%lg  %lg\n", p, mlik_interp(mlik, p) );
    }
  fclose(fp);
  #endif

  return TRUE;
  }



double mlik_interp(Interpar *mlik, double p){
  double m, dp;

  if(p<=mlik->min)  {
    dp=p-mlik->min; 
    m=myinterp(mlik, mlik->min)+dp*mlik->slop_min;
    }
  else if(p>=mlik->max)  {
    dp=p-mlik->max; 
    m=myinterp(mlik, mlik->max)+dp*mlik->slop_max;
    }
  else{
    m=myinterp(mlik, p); 
    }

  return m;
  }

/*->> end of max likelihood curve subroutines <<-*/



void load_stat_disp(float *disp, float *disp_phi, float *disp_model, float *phi, 
                    float *phi_model, char *fname, long long npart) {
  FILE *fp;
  if(!(fp=fopen(fname, "r"))) {
    printf("can't open file `%s`\n", fname);
    exit(0);
    }

  //->>  loading <<- //
  float *tmp;
  tmp=(float *)fftwf_malloc(sizeof(float)*npart*3);

  //->>
  fread(disp, sizeof(float), npart*3, fp);
  fread(disp_model, sizeof(float), npart*3, fp);
  fread(tmp, sizeof(float), npart*3, fp);

  fread(tmp, sizeof(float), npart, fp);
  fread(phi, sizeof(float), npart, fp);
  fread(disp_phi, sizeof(float), npart*3, fp);
  fread(tmp, sizeof(float), npart, fp);
  fread(phi_model, sizeof(float), npart, fp);

  //->> 
  fclose(fp);
  fftwf_free(tmp);
  return;
  }




void phi_maximum_fitting(SimInfo *s, Interpar *mlik, float *phi_model, 
                        float *phi_nl, float *phi_cb, long long ngrid)  {
  //->> return the maximum of <<- //
  long long i, j, k;
  float p;

  #ifdef _OMP_
  #pragma omp parallel for private(i,j,k,p)
  #endif
  for(i=0; i<ngrid; i++)
    for(j=0; j<ngrid; j++)
      for(k=0; k<ngrid; k++) {
        p=ArrayAccess3D(phi_model, ngrid, i, j, k);
        ArrayAccess3D(phi_nl, ngrid, i, j, k)=mlik_interp(mlik, p);
        ArrayAccess3D(phi_cb, ngrid, i, j, k)=p+ArrayAccess3D(phi_nl, ngrid, i, j, k);
        }

  return;
  }






void phi_mlik_displacement(SimInfo *s, Pdata_pos *p, Interpar *mlik, float *disp, 
          float *disp_phi, float *disp_model, float *phi, float *phi_model, 
          long long ngrid, double boxsize, char *stat_disp_fname, int import_disp) {
  /* ->> wrapper for maximum likelihood reconstruction <<- */

  long long i, j, k, npart;
  npart=ngrid*ngrid*ngrid;

  // ->> load displacement field first <<- //
  if(import_disp==TRUE) {
    disp=(float *)fftwf_malloc(sizeof(float)*npart*3);
    disp_phi=(float *)fftwf_malloc(sizeof(float)*npart*3);
    disp_model=(float *)fftwf_malloc(sizeof(float)*npart*3);
    phi=(float *)fftwf_malloc(sizeof(float)*npart);
    phi_model=(float *)fftwf_malloc(sizeof(float)*npart);

    load_stat_disp(disp, disp_phi, disp_model, phi, phi_model, stat_disp_fname, npart);
    }

  // ->> maximum of phi fitting <<- //
  float *phi_nl, *phi_cb, *disp_rec, *d_rec; 
  phi_nl=(float *)fftwf_malloc(sizeof(float)*npart);
  phi_cb=(float *)fftwf_malloc(sizeof(float)*npart);

  disp_rec=(float *)fftwf_malloc(sizeof(float)*npart*3);
  d_rec=(float *)fftwf_malloc(sizeof(float)*npart);

  phi_maximum_fitting(mlik, phi_model, phi_nl, phi_cb, ngrid);

  // ->> obtain the reconstructed displacement field from phi <<- //
  fft_gradient(phi_cb, disp_rec, boxsize, ngrid);

  //->> move particles <<- //
  long long ngrid_xyz[3];
  double pmin[3], pmax[3], dpart[3];

  for(i=0; i<3; i++) {ngrid_xyz[i]=ngrid;}
  get_particle_boundary(p, boxsize, npart, ngrid_xyz, pmin, pmax, dpart);

   
  p_disp=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
  general_particle_mover(p, p_disp, disp_rec, boxsize, ngrid, FALSE); //->> do not interpolate <<- //
 
  //->> reconstructed density field <<- //
  
   //->> 
   printf("chnage move_grid_general");
   abort();
   move_grid_general(s, , disp_rec, , si, boxsize, ngrid, FALSE);
  


  //->> output reconstructed displacement & density <<- //
  output_maxlikelihood_data(float *phi, float *phi_model, float *phi_nl);


  // ->> free <<- // 
  fftwf_free(disp); fftwf_free(disp_model); fftwf_free(disp_phi);
  fftwf_free(phi); fftwf_free(phi_model);  
  fftwf_free(phi_nl); fftwf_free(phi_cb); fftwf_free(d_rec);

  return;
  }











void output_maxlikelihood_data(SimInfo *s, float *disp, float *disp_rec, float *d_rec, long long npart){

  return;
  }



