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


#ifdef _MPI_
  #include "mpi.h"
#endif

#ifdef _OMP_
  #include <omp.h>
#endif



void load_stat_disp(float *disp, float *disp_model, float *phi, 
                                float *phi_model, char *fname){
  FILE *fp;
  if(!(fp=fopen(fname, "r"))) {
    printf("can't open file `%s`\n", fname);
    exit(0);
    }
  //->>  

  float *tmp;

  fread(disp);

  return;
  }




void phi_mlik_displacement(SimInfo *s, Interpar *mlik, float *disp, float *phi,
                           float *disp_model, float *phi_model, long long ngrid, 
	                   char *stat_disp_fname, int import_disp) {
  long long i, j;

  // ->> load displacement field first <<- //
  load_stat_disp(disp, disp_model, phi, phi_model, stat_disp_fname);

  //->>  




  //->> output reconstructed displacement <<- //
   

  return;
  }








void output_maxlikelihood_data(){

  return;
  }
