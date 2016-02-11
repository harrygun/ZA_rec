#ifndef _H_HISTOGRAM_
#define _H_HISTOGRAM_

  #include <gsl/gsl_histogram2d.h>

  typedef struct {
    gsl_histogram2d *h;
    gsl_histogram2d_pdf *p; 
  
    int grid[2];
    double boundary[2][2];
    } Histopar2d;


#endif
