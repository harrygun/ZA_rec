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



void transfer_func_init( ) {
  // ->> import transfer function and interpolate <<- //


  return;
  }





void get_stat_disp_model(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *stat_disp_model_type) {
  // ->>  <<- //
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
