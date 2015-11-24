  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include<fftw3.h>


  #define pi M_PI



void poisson_solver(double ***d, double ***phi, int ngrid)  {
  /* ->> Poisson Solver with FFT <<- */

  int l, m, n;
  double kx, ky, kz, sin2x, sin2y, sin2z, greens;
  double fac=1.0/(double)(ngrid*ngrid*ngrid);

  fftw_complex dk[ngrid][ngrid][ngrid/2+1];
  fftw_plan pforward, pbackward;
  
  pforward = fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid, &d[0][0][0], &dk[0][0][0], FFTW_ESTIMATE);
  fftw_execute(pforward);  
  
  /* work out the green's function and the FFT of phi*/
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid/2+1; n++){
        kx = 2.*pi*l/(double)ngrid;
        ky = 2.*pi*m/(double)ngrid;
        kz = 2.*pi*n/(double)ngrid;
        
        sin2x = sin(kx/2.)*sin(kx/2.);
        sin2y = sin(ky/2.)*sin(ky/2.);
        sin2z = sin(kz/2.)*sin(kz/2.);
        
        if ((l==0) && (m==0) && (n==0)) greens = 0.;
        else greens = - (double)3./8./a/(sin2x+sin2y+sin2z);
      	    
        dk[l][m][n][0] *= greens;
        dk[l][m][n][1] *= greens;
        }
  
  /* find the inverse FFT of phi */
  pbackward = fftw_plan_dft_c2r_3d(ngrid, ngrid, ngrid, &dk[0][0][0], &phi[0][0][0], FFTW_ESTIMATE);
  fftw_execute(pbackward);
  
  for (l=0; l<ngrid; l++)
    for (m=0; m<ngrid; m++)
      for (n=0; n<ngrid; n++)
      	phi[l][m][n] *= fac;
  
  fftw_destroy_plan(pforward);
  fftw_destroy_plan(pbackward);

  return;
  }
