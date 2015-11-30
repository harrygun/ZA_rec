  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"



void move_particle(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, float *s, int s_intp){
  //-> move particles <<- //
  int i, j;

  if(s_intp!=TRUE){
    // ->> do not interpolate <<- //

    for(i=0; i<s->npart; i++) {
       
      }

    }

  else{ }

  return;
  }
