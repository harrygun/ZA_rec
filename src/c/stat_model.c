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
  get_real_displacement(s, pinit, pinit, disp_init, disp_calmethod);

  // ->>  get final displacement <<- //
  get_real_displacement(s, p, pinit, disp, disp_calmethod);

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



Interpar *transfer_func_init(char *fname) {
  // ->> import transfer function and interpolate <<- //
  Interpar *tf= (Interpar *)malloc(3*sizeof(Interpar));
  FILE *fp = fopen(fname, "r"); 

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

    //printf("tf boundary:  %lg   %lg  %lg  %lg\n", tf[j].min, tf[j].max, tf[j].slop_min, tf[j].slop_max);
    }

  fclose(fp);


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

  return tf;
  }



void transfer_func_finalize(Interpar *tf){
  int i;

  for(i=0; i<3; i++)
    myinterp_free(&tf[i]);

  return;
  }



/* -> displacement field manipulation  <- */
void load_displacement(Cospar *cp, SimInfo *s, Pdata_pos *p, Pdata_pos *pinit, 
                            float *disp, float *disp_lpt)  {

  // -> get the stochastic term of displacement <<- //
  int i;
  char *disp_calmethod="grid_wise";


  // ->> import initial displacement <<- //
  Pdata_pos *pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata_pos));
  load_cita_simulation_position(fname_part_init, pinit, s->npart);
 
  // ->>  get initial displacement <<- //
  get_real_displacement(s, pinit, pinit, disp_init, disp_calmethod);

  // ->>  get final displacement <<- //
  get_real_displacement(s, p, pinit, disp, disp_calmethod);

  free(pinit);
  return;
  }



void get_disp_mc(Cospar *cp, SimInfo *s, Interpar *tf){
  // ->> free <<- //
  return;
  }



void get_stat_disp_model(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *stat_disp_model_type) {
  // ->> obtain statistical displacement model <<- //
  float *disp, *disp_model;
  Pdata_pos *pinit;

  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  disp_model=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  // ->> obtain real displacement <<- //
  char *disp_calmethod="grid_wise";

  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
    // ->> load initial position <<- //
    pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
    load_cita_simulation_position(fname_part_init, pinit, s->npart);
    }
  get_real_displacement(s, p, pinit, disp, disp_calmethod);


  // ->> obtain model displacement <<- //
  char *model_type;
  if(stat_disp_model_type==NULL) {
    sprintf(model_type, "ZA");
    }
  else {
    model_type=stat_disp_model_type;
    }
  get_model_displacement(s, p, d, disp_model, model_type);


  // ->> construct model <<- //





  free(disp);
  return;
  }
