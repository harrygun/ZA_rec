  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"



void move_particle(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, float *s, int s_intp){
  //-> move particles <<- //
  int i, j, ip;


  for(ip=0; ip<s->npart; ip++) {

    if(s_intp!=TRUE){
      // ->> do not interpolate <<- //
      }
    else{  // neglect position difference
      for(i=0; i<3; i++)
        p[ip].pos[i]   

      }

    }


  return;
  }
