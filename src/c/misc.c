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

  #ifdef _OMP_
  #include <omp.h>
  #endif


void mat_inv_3d_float(float mat[3][3], float imat[3][3]){

  //->> 
  imat[0][0]=(-mat[1][2]*mat[2][1]+mat[1][1]*mat[2][2])/(-mat[0][2]*mat[1][1]*mat[2][0] 
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[0][1]=(mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[0][2]=(-mat[0][2]*mat[1][1] + mat[0][1]*mat[1][2])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);

  //->> 
  imat[1][0]=(mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[1][1]=(-mat[0][2]*mat[2][0] + mat[0][0]*mat[2][2])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[1][2]=(mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
 
  //->> 
  imat[2][0]=(-mat[1][1]*mat[2][0] + mat[1][0]*mat[2][1])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[2][1]=(mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);
  imat[2][2]=(-mat[0][1]*mat[1][0] + mat[0][0]*mat[1][1])/(-mat[0][2]*mat[1][1]*mat[2][0]
             + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] 
	     - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] 
	     + mat[0][0]*mat[1][1]*mat[2][2]);

  return;
  }




void mat_multiply_3d_float(float mat[3][3], float x[3], float y[3]){
  int i, j;

  for(i=0; i<3; i++) {
    y[i]=0.;
    for(j=0; j<3; j++)
      y[i]+=mat[i][j]*x[j];
    }
  return;
  }


// ->> trim the boundary of cubic <<- //
void cubic_trim(float *mi, float *mo, size_t ngrid, size_t ntrim){
  size_t i, j, k, nt;

  nt=ngrid-2*ntrim; 

  #ifdef _OMP_
  #pragma omp parallel for private(i,j,k)
  #endif
  for(i=0; i<nt; i++)
    for(j=0; j<nt; j++)
      for(k=0; k<nt; k++) {
        mo[] = mi ;

        }

  return;
  }

