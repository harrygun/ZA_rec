  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"



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
      prinf("particle-moving interpolation NOT supported yet.\n");
      abort();
      }
    else{ abort(); }

    }

  return;
  }
