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







//->> initialization of 2d histogram  <<- //
void histogram2d_init(Histopar2d *his) {
  // ->> initialization of histogram <<- //

  his->h=gsl_histogram2d_alloc(his->grid[0], his->grid[1]);
  his->p=gsl_histogram2d_pdf_alloc(his->h->nx, his->h->ny);

  // 
  gsl_histogram2d_pdf_init(his->p, his->h);

  return;
  }



void histogram2d_free(Histopar2d *his) {

  gsl_histogram2d_pdf_free(his->p);
  gsl_histogram2d_free(his->h);
  return;
  }



//void 
