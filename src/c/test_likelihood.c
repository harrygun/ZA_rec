  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
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
  #include "poisson.h"
  #include "reconstruction_partmoving.h"
  #include "fourier.h"

  #include "stat_model.h"



#ifdef _OMP_
  #include <omp.h>
#endif




void test_displacement(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                       char *fname_out) {
  // ->> try to build the statistical model from real displacement <<- //
  //get_stat_disp_model(s, p, d, fname_part_init, NULL);
  int i;
  float *disp, *disp_model;

  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  disp_model=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  // ->> obtain real displacement <<- //
  get_real_displacement(s, p, disp, fname_part_init);

  // ->> smooth the field <<- //
  for(i=0; i<3; i++)
    smooth_field(&disp[i*s->ngrid*s->ngrid*s->ngrid], s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R);

  // ->> obtain model displacement <<- //
  char *model_type="ZA"; //"2LPT";  //"ZA";
  get_model_displacement(s, p, d, disp_model, model_type);

  // ->> construct model <<- //
  FILE *fp=fopen(fname_out, "wb");

  fwrite(disp, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(disp_model, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);

  fclose(fp);

  free(disp); free(disp_model);
  return;
  }


