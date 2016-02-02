
  #ifndef _H_MY_INTERPOLATE_
     #define _H_MY_INTERPOLATE_

      #include <gsl/gsl_spline.h>
      #include <gsl/gsl_integration.h>
  
      typedef struct {

        gsl_interp_accel *acc;
        gsl_spline *spl;
  
        } Interpar;
  


      int myinterp_init( Interpar *f, double *x, double *y, int n);
      //Interpar * gslinterp_init( double *x, double *y, int n) ;
      int gslinterp_init(Interpar *f, double *x, double *y, int n);
      void * gslinterp_void_init(double *x, double *y, int n) ;


      int myinterp_free( Interpar *f);

      double myinterp(Interpar *f, double x);



      double bilinear_interp_ele(double f11, double f12, double f21, double f22, 
             double x1, double x2, double y1, double y2, double x, double y);


  #endif
