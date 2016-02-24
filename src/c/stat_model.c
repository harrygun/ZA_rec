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

  #include "displacement.h"


  #ifdef _OMP_
  #include <omp.h>
  #endif




void disp_field_tranfunc_precal(SimInfo *s, Pdata_pos *p, float *d, 
                                 char *fname_part_init, char *fname_out) {
  // ->> output displacement field for calculating transfer function <<- //
  int i;
  float *disp_init, *disp;
  disp_init=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  // ->> get displacement field differently <<- //
  Pdata_pos *pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata_pos));
  load_cita_simulation_position(fname_part_init, pinit, s->npart);

  // ->>  get initial displacement <<- //
  char *disp_calmethod="grid_wise";
  get_real_displacement(s, pinit, pinit, disp_init, s->drift_init, disp_calmethod, 1.);

  // ->>  get final displacement <<- //
  get_real_displacement(s, p, pinit, disp, s->drift, disp_calmethod, 1.);

  // ->> output <<- //
  FILE *fp=fopen(fname_out, "wb");
  fwrite(disp_init, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(disp, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fclose(fp);

  // ->> free <<- //
  free(disp_init); free(disp);
  free(pinit);
  return;
  }



double tk_interp(Interpar *tf, double k){

  double tk_b, tk, dk;

  if(k<=tf->min)  {
    tk_b=myinterp(tf, tf->min);
    dk=k-tf->min; 
    tk=exp(log(tk_b)+dk*tf->slop_min);
    }

  else if(k>=tf->max){
    tk_b=myinterp(tf, tf->max);
    dk=k-tf->max; 
    tk=exp(log(tk_b)+dk*tf->slop_max);
    }

  else{
    tk=myinterp(tf, k); 
    }

  if(tk<0)  tk=0.;
  if(tk>1)  tk=1.;

  return tk;
  }



//Interpar *transfer_func_init(char *fname) {
int transfer_func_init(Interpar *tf, char *fname) {
  // ->> import transfer function and interpolate <<- //
  FILE *fp;

  if(!(fp=fopen(fname, "r"))) {
    return FALSE; 
    }

  // ->> importing data <<- //
  int i, j, line, extdidx=15;
  line=countFILEline(fp);

  double *k, *tfk, kk, dtf, dk;
  k=(double *)malloc(3*line*sizeof(double));
  tfk=(double *)malloc(3*line*sizeof(double));

  for(i=0; i<line; i++) {
    for(j=0; j<3; j++) {
      fscanf(fp, "%lg  %lg ", &ArrayAccess2D_n2(k, 3, line, j, i), &ArrayAccess2D_n2(tfk, 3, line, j, i));
      }
    fscanf(fp, "\n");
    }

  // ->> initialize interpolator <<- //
  for(j=0; j<3; j++) {
    myinterp_init(&tf[j], &ArrayAccess2D_n2(k, 3, line, j, 0), &ArrayAccess2D_n2(tfk, 3, line, j, 0), line);

    //->> extrapolation initialization <<- //
    sprintf(tf[j].extra_type, "log-linear");

    tf[j].min=ArrayAccess2D_n2(k, 3, line, j, 0);
    tf[j].max=ArrayAccess2D_n2(k, 3, line, j, line-1);

    dk=ArrayAccess2D_n2(k, 3, line, j, extdidx)-ArrayAccess2D_n2(k, 3, line, j, 0);
    dtf=log(ArrayAccess2D_n2(tfk, 3, line, j, extdidx))-log(ArrayAccess2D_n2(tfk, 3, line, j, 0));
    tf[j].slop_min=dtf/dk;

    dk=ArrayAccess2D_n2(k, 3, line, j, line-extdidx/2)-ArrayAccess2D_n2(k, 3, line, j, line-extdidx);
    dtf=log(ArrayAccess2D_n2(tfk, 3, line, j, line-extdidx/2))-log(ArrayAccess2D_n2(tfk, 3, line, j, line-extdidx));
    tf[j].slop_max=dtf/dk;

    printf("tf boundary [%d]:  %lg   %lg  %lg  %lg\n", j, tf[j].min, tf[j].max, tf[j].slop_min, tf[j].slop_max);
    }

  fclose(fp);


  //#define _DO_TF_TEST_
  #ifdef _DO_TF_TEST_
  fp=fopen("result/test_tf.dat", "w");
  for(i=0; i<line; i++) {
    for(j=0; j<3; j++) {
      //kk=ArrayAccess2D_n2(k, 3, line, j, i);
      kk=pow(10., -2+(double)i*4./((double)line) ); 
      printf("%lg  %lg  ", kk, tk_interp(&tf[j], kk) );
      fprintf(fp, "%lg  %lg  ", kk, tk_interp(&tf[j], kk) );
      //printf("%lg  %lg  ", ArrayAccess2D_n2(k, 3, line, j, i), ArrayAccess2D_n2(tfk, 3, line, j, i));
      }
    printf("\n");
    fprintf(fp, "\n");
    } 
  fclose(fp);
  #endif


  free(k); free(tfk);

  return TRUE;
  }



// ->> free interpolator <<- //
void transfer_func_finalize(Interpar *tf){
  int i;
  for(i=0; i<3; i++)
    myinterp_free(&tf[i]);

  return;
  }



/* ->> max likelihood curve <<- */
int phi_mlik_init(Interpar *mlik, char *fname){

  FILE *fp;
  if(!(fp=fopen(fname, "r"))) {
    return FALSE; }

  // ->> importing data <<- //
  int i, j, line;
  line=countFILEline(fp);

  double *phi, *mlik_phi, *smlik_phi, p, ml;
  phi=(double *)malloc(line*sizeof(double));
  mlik_phi=(double *)malloc(line*sizeof(double));
  smlik_phi=(double *)malloc(line*sizeof(double));

  for(i=0; i<line; i++) 
    fscanf(fp, "%lg  %lg  %lg\n", &phi[i], &mlik_phi[i], &smlik_phi[i]);


  // ->> initialize interpolator <<- //
  myinterp_init(mlik, phi, smlik_phi, line);

  mlik->min=phi[0];
  mlik->max=phi[line-1];



  #define _DO_MLIK_TEST_
  #ifdef _DO_MLIK_TEST_
  fp=fopen("result/test_mlik.dat", "w");

  for(i=0; i<500; i++) {
    p=-300+i*600./500.; 
    printf("%lg  %lg\n", p, mlik_interp(mlik, p) );
    fprintf(fp, "%lg  %lg\n", p, mlik_interp(mlik, p) );
    }
  fclose(fp);
  #endif

  abort();


  return TRUE;
  }



double mlik_interp(Interpar *mlik, double p){
  double p, m, dp;

  if(p<=mlik->min)  {
    dp=p-mlik->min; 
    m=myinterp(mlik, mlik->min)+dp*mlik->slop_min;
    }

  else if(k>=tf->max){
    tk_b=myinterp(tf, tf->max);
    dk=k-tf->max; 
    m=exp(log(tk_b)+dk*tf->slop_max);
    }

  else{
    m=myinterp(mlik, p); 
    }

  return m;
  }




/* -> displacement field manipulation  <- */
void load_displacement(Cospar *cp, SimInfo *s, Pdata_pos *p, float *disp, 
               float *disp_lpt, char *fname_part_init, char *fname_pid_init)  {
  // -> get the stochastic term of displacement <<- //
  double fac;
  //char *disp_calmethod="grid_PID";  //"grid_wise";
  char *disp_calmethod="init_pos_PID";  //"grid_wise";

  // ->> rescaling factor for initial data <<- //
  fac=Dp(cp, cp->z)/Dp(cp, cp->zinit);

  // ->> import initial displacement <<- //
  Pdata_pos *pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata_pos));
  load_cita_simulation_position_pid(fname_part_init, fname_pid_init, pinit, s->npart);

  // ->>  get final displacement <<- //
  get_real_displacement(s, p, pinit, disp, s->drift, disp_calmethod, 1.);
 
  // ->>  get initial displacement <<- //
  disp_calmethod="grid_wise";
  get_real_displacement(s, pinit, pinit, disp_lpt, s->drift_init, disp_calmethod, fac);

  free(pinit);
  return;
  }




void disp_stat_separation(Cospar *cp, float *disp, float *disp_lpt, float *disp_mc, 
                      Interpar *tf, float boxsize, long long npart, long long ngrid)  {
  // ->> For each particle, statistical separate the into  <<- //
  //     deterministic and stochastic contributions .      <<- //
  long long ip, i, j, k;
  
  //->> convolve displacement field with transfer function <<- //
  for(i=0; i<3; i++) {
    smooth_field(&ArrayAccess2D_n2(disp_lpt, 3, npart, i, 0), boxsize, ngrid, 
                            _ANISOTROPIC_INTERPOLATOR_SMOOTH_, 0., &tf[i]);
    }

  // ->> now get the mode-coupling term <<- //
  #ifdef _OMP_
  #pragma omp parallel for private(ip,i)
  #endif
  for(ip=0; ip<npart; ip++) 
    for(i=0; i<3; i++)    {
      ArrayAccess2D_n2(disp_mc, 3, npart, i, ip)=ArrayAccess2D_n2(disp, 3, npart, i, ip)
                                                   -ArrayAccess2D_n2(disp_lpt, 3, npart, i, ip);
    }

  return;
  }





void get_stat_disp_model(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *stat_disp_model_type) {
  // ->> construct model <<- //

  return;
  }




/* ->> output fields <<- */
void output_real_disp_field(float *disp, float *disp_lpt, size_t ngrid, char *fname_out){
  //->> only write full and LPT displacement field <<- //

  printf("output displacement field:  %s\n", fname_out); 
  fflush(stdout);

  FILE *fp=fopen(fname_out, "wb");

  fwrite(disp, sizeof(float), ngrid*ngrid*ngrid*3, fp);
  fwrite(disp_lpt, sizeof(float), ngrid*ngrid*ngrid*3, fp);

  fclose(fp);
  return;
  }



/* ->> output routines <<- */
void output_stat_disp_model(float *disp, float *disp_lpt, float *disp_mc, 
                      size_t ngrid, size_t ngrid_trimmed, char *fname_out){
  // ->> write displacement field into files <<- //
  //
  FILE *fp=fopen(fname_out, "wb");
  fwrite(disp, sizeof(float), ngrid*ngrid*ngrid*3, fp);
  fwrite(disp_lpt, sizeof(float), ngrid*ngrid*ngrid*3, fp);
  fwrite(disp_mc, sizeof(float), ngrid*ngrid*ngrid*3, fp);

  fclose(fp);
  return;
  }




void output_stat_disp_potential_model(float *disp, float *disp_lpt, float *disp_mc, 
                float *div, float *phi, float *disp_phi, float *div_lpt, float *phi_lpt, 
                float *disp_phi_lpt,  size_t ngrid, size_t ngrid_trim, char *fname_out){
  // ->> write displacement field into files <<- //
  //
  size_t ng, tng;
  ng=ngrid*ngrid*ngrid;
  tng=ngrid_trim*ngrid_trim*ngrid_trim;

  FILE *fp=fopen(fname_out, "wb");

  //fwrite(disp, sizeof(float), ng*3, fp);
  //fwrite(disp_lpt, sizeof(float), ng*3, fp);
  //fwrite(disp_mc, sizeof(float), ng*3, fp);
  //
  fwrite(disp, sizeof(float), tng*3, fp);
  fwrite(disp_lpt, sizeof(float), tng*3, fp);
  fwrite(disp_mc, sizeof(float), tng*3, fp);

  fwrite(div, sizeof(float), tng, fp);
  fwrite(phi, sizeof(float), tng, fp);
  fwrite(disp_phi, sizeof(float), tng*3, fp);
  fwrite(div_lpt, sizeof(float), tng, fp);
  fwrite(phi_lpt, sizeof(float), tng, fp);

  if(disp_phi_lpt!=NULL)
    fwrite(disp_phi_lpt, sizeof(float), tng*3, fp);


  fclose(fp);
  return;
  }
