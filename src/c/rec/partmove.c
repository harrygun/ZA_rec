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
  float moved_pos, xmin, xmax, dx;

  // ->> box boundary <<- //
  xmin=0.; xmax=s->boxsize;
  dx=(xmax-xmin)/(float)s->ngrid;
       
  printf("->> displaceing particles...\n");

  for(ip=0; ip<s->npart; ip++) {

    // ->> do not interpolate <<- //
    if(s_intp==FALSE){
      for(i=0; i<3; i++){
        moved_pos=p[ip].pos[i]+ArrayAccess2D_n2(si, 3, s->npart, i, ip);
        // ->> periodic boundary condition <<- //
	if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
	else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
	else{
	  moved[ip].pos[i]=moved_pos; }
	}
      }

    // ->> interpolate shift field onto particle position <<- //
    else if(s_intp==TRUE){  
      printf("particle-moving interpolation NOT supported yet.\n");
      abort();
      }
    else{abort();}

    }

  printf("->> particles displacement is done.\n");
  return;
  }



void move_grid(SimInfo *s, Pdata_pos *moved, float *si, int s_intp){
  //->> grid moving, no need to generate the grid <<- //
  int i, j, k, m, ip;
  float grid[3], xmin, xmax, dx, moved_pos;
       
  // ->> box boundary <<- //
  xmin=0.; xmax=s->boxsize;
  dx=(xmax-xmin)/(float)s->ngrid;

  printf("\n->> shifting uniform grid ...\n");

  for(i=0; i<s->ngrid; i++)
    for(j=0; j<s->ngrid; j++)
      for(k=0; k<s->ngrid; k++) {

        // ->> do not interpolate <<- //
        if(s_intp==FALSE){
          // ->> grid index <<- //
          ip=MemIdx3D(s->ngrid, i, j, k);

	  grid[0]=xmin+i*dx;
	  grid[1]=xmin+j*dx;
	  grid[2]=xmin+k*dx;

	  for(m=0; m<3; m++)  {
            moved_pos=grid[m]+ArrayAccess2D_n2(si, 3, s->npart, m, ip);
            // ->> periodic boundary condition <<- //
	    if(moved_pos<0){
              moved[ip].pos[m]=moved_pos+xmax; }
	    else if(moved_pos>=xmax){
              moved[ip].pos[m]=moved_pos-xmax; }
	    else{
	      moved[ip].pos[m]=moved_pos; }
            }

          }

        // ->> interpolate shift field onto particle position <<- //
        else if(s_intp==TRUE){  
          printf("grid-moving interpolation NOT supported yet.\n");
          abort();
          }
        else{abort();}
      }

  printf("->> shifting grid is done.\n");
  return;
  }
