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

  #include "const.h"
  #include "parvar.h"
  #include "io.h"
  #include "cic.h"
  #include "za_reconstruction.h"



void move_particle(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, float *si, int s_intp){
  //-> move particles <<- //
  int i, ip;

       
  for(ip=0; ip<s->npart; ip++) {

    if(s_intp==FALSE){
      // ->> do not interpolate <<- //
      for(i=0; i<3; i++) {
        moved[ip].pos[i]=p[ip].pos[i]+ArrayAccess2D_n2(si, 3, s->npart, i, ip);
	}
      }
    else if(s_intp==TRUE){  
      // ->> interpolate shift field onto particle position <<- //
      printf("particle-moving interpolation NOT supported yet.\n");
      abort();
      }
    else{abort();}

    }

  return;
  }



void move_grid(SimInfo *s, Pdata_pos *moved, float *si, int s_intp){
  //->> grid moving, no need to generate the grid <<- //
  int i, ip;
  float grid[3];
       
  for(ip=0; ip<s->npart; ip++) {

    if(s_intp==FALSE){
      // ->> do not interpolate <<- //
      for(i=0; i<3; i++) {
        grid[i]= ;
        moved[ip].pos[i]=grid[i]+ArrayAccess2D_n2(si, 3, s->npart, i, ip);
	}
      }
    else if(s_intp==TRUE){  
      // ->> interpolate shift field onto particle position <<- //
      printf("grid-moving interpolation NOT supported yet.\n");
      abort();
      }
    else{abort();}

    }


  return;
  }
