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



void get_real_displacement(SimInfo *s, Pdata_pos *p, float *disp, 
                           char *fname_part_init, char *disp_calmethod) {
  // ->> get real displacement field from simulation directly <<- //

  long long ip, i;
  Pdata_pos *pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
  load_cita_simulation_position(fname_part_init, pinit, s->npart);


  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
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

    // ->> re-arrange data <<- //
    #ifdef _OMP_
    #pragma omp parallel for private(i,j,k,m,ip,grid,moved_pos)
    #endif
    for(i=0; i<s->ngrid; i++)
      for(j=0; j<s->ngrid; j++)
        for(k=0; k<s->ngrid; k++) {

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

  }



  free(pinit);
  return;
  }





void get_model_displacement(SimInfo *s, Pdata_pos *p, float *d, float *disp, char *model_disp_type){
    

  if(strcmp(model_disp_type, "ZA")==0 ) {
    za_displacement(s, d, disp);
    }
  else if(strcmp(model_disp_type, "2LPT")==0 ) {
    displacement_2lpt(s, d, disp);
    }
  else abort();

  return;
  }



