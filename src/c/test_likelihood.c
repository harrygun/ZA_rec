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

  #include "stat_model.h"



#ifdef _OMP_
  #include <omp.h>
#endif




void test_displacement(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init) {


  // ->> try to build the statistical model from real displacement <<- //
  get_stat_disp_model(s, p, d, fname_part_init, NULL);

  // ->>  <<- //




  return;
  }
