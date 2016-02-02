#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "interpolation.h"


  int myinterp_init( Interpar *f, double *x, double *y, int n)  {


      f->acc = gsl_interp_accel_alloc();
      f->spl = gsl_spline_alloc(gsl_interp_cspline, n);
 
      gsl_spline_init(f->spl, x, y, n);

      return TRUE;
    }



 int gslinterp_init(Interpar *f, double *x, double *y, int n)  {

      f=(Interpar *)malloc(sizeof(Interpar));

      f->acc = gsl_interp_accel_alloc();
      f->spl = gsl_spline_alloc(gsl_interp_cspline, n);
 
      gsl_spline_init(f->spl, x, y, n);

      printf("gsl_interpolater initialization done.\n"); fflush(stdout);
      printf("n=%d\n", n);
      double xx, yy;
      xx=1.0;
      yy=myinterp(f, xx);

      printf("INSIDE:  x=%lg, y=%lg\n", xx, yy);

      return TRUE;
    }


 void * gslinterp_void_init(double *x, double *y, int n)  {

      Interpar *f=(Interpar *)malloc(sizeof(Interpar));

      f->acc = gsl_interp_accel_alloc();
      f->spl = gsl_spline_alloc(gsl_interp_cspline, n);
 
      gsl_spline_init(f->spl, x, y, n);

      return (void *)f;
    }




  int myinterp_free( Interpar *f)  {

      gsl_spline_free(f->spl);
      gsl_interp_accel_free(f->acc);

      return TRUE;
    }


  double myinterp(Interpar *f, double x)  {

       //printf("here is ok, x=%lg\n", x); fflush(stdout);
       double y= gsl_spline_eval(f->spl, x, f->acc);
       //printf("here is ok, y=%lg\n", y); fflush(stdout);

       return y;
    }



  double bilinear_interp_ele(double f11, double f12, double f21, double f22, 
             double x1, double x2, double y1, double y2, double x, double y) {
    double fxy;

    fxy = f11/(x2-x1)/(y2-y1)*(x2-x)*(y2-y) +
          f21/(x2-x1)/(y2-y1)*(x-x1)*(y2-y) +
          f12/(x2-x1)/(y2-y1)*(x2-x)*(y-y1) +
          f22/(x2-x1)/(y2-y1)*(x-x1)*(y-y1);
    return fxy;
  }


