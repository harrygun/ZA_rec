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

  #include "parvar.h"
  #include "io.h"
  #include "cic.h"
  #include "misc.h"
  #include "poisson.h"
  #include "fourier.h"

  #include "backward_displacement.h"


  #ifdef _OMP_
  #include <omp.h>
  #endif



void get_real_displacement(SimInfo *s, Pdata_pos *p, Pdata_pos *pinit, float *disp, 
                           char *disp_calmethod) {
  // ->> get real displacement field from simulation directly <<- //

  long long ip, i, j, k, m;
  float grid[3], xmin, xmax, dx;

  // ->> box boundary <<- //
  xmin=0.; xmax=s->boxsize;
  dx=(xmax-xmin)/(float)s->ngrid;

  // ->>
  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
    printf("direct subtraction for real_displacement.\n"); fflush(stdout);

    #ifdef _OMP_
    #pragma omp parallel for private(ip,i)
    #endif
    for(ip=0; ip<s->npart; ip++) {
      for(i=0; i<3; i++) {
        ArrayAccess2D_n2(disp, 3, s->npart, i, ip)=p[ip].pos[i]-pinit[ip].pos[i];
        }
      }
    }
  else if( strcmp(disp_calmethod, "grid_wise")==0 ){
    printf("grid-wise calculation for real_displacement.\n"); fflush(stdout);

    // ->> re-arrange data into grid <<- //
    #ifdef _OMP_
    #pragma omp parallel for private(i,j,k,m,ip,grid)
    #endif
    for(i=0; i<s->ngrid; i++)
      for(j=0; j<s->ngrid; j++)
        for(k=0; k<s->ngrid; k++) {

          // ->> grid index <<- //
          ip=MemIdx3D(s->ngrid, i, j, k);

          grid[0]=xmin+i*dx;
          grid[1]=xmin+j*dx;
          grid[2]=xmin+k*dx;

          for(m=0; m<3; m++){
            //ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=p[ip].pos[m]-pinit[ip].pos[m];
            //ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=p[ip].pos[2-m]-pinit[ip].pos[2-m];
	    
            ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=p[ip].pos[m]-grid[m];
            //ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=p[ip].pos[2-m]-grid[m];
            }

          }
    }


  return;
  }





void get_model_displacement(SimInfo *s, Pdata_pos *p, float *d, float *disp, char *fname_part_init, char *model_disp_type){

  int model_from_init_pos=TRUE;

  Pdata_pos *pinit;
  double dmean;

  if (model_from_init_pos==TRUE){
    pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
    load_cita_simulation_position(fname_part_init, pinit, s->npart);
      
    // ->> estimate the density field from initial field <<- //
    dmean=cic_density(pinit, d, s->boxsize, s->particle_mass, s->npart, s->ngrid_xyz, s); 
    }


  if(strcmp(model_disp_type, "ZA")==0 ) {
    za_displacement(s, d, disp);
    }
  else if(strcmp(model_disp_type, "2LPT")==0 ) {
    displacement_2lpt(s, d, disp);
    }
  else abort();


  if (model_from_init_pos==TRUE){
    free(pinit); }

  return;
  }



