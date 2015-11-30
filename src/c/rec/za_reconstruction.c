  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "parvar.h"
  #include "poisson.h"
  #include "cic.h"
  #include "partmove.h"



void za_displacement(SimInfo *s, float *d, float *disp) {
  // ->> obtain ZA displacement field from density field <<- // 

  int fft_return_type;
  float *phi, phi_ij;

  fft_return_type=_RETURN_GRADIENT_;
  phi=(float *)fftwf_malloc(sizeof(float)*s->ngrid_xyz[0]*s->ngrid_xyz[1]*s->ngrid_xyz[2]);

  // ->> solve FFTW  <<- //
  printf("\n->> Solve Poisson equation with FFT.\n");
  poisson_solver_float(d, phi, disp, phi_ij, s->boxsize, s->ngrid, 
                       s->smooth_type_flag, s->smooth_R, fft_return_type);
  printf("->> FFT is Done.\n");

  fftwf_free(phi);
  return;
  }





void za_reconstruction(SimInfo *s, Pdata_pos *p, float *d, float *drec, char *rec_type){
  // ->> Performing ZA reconstruction  <<- //
  int disp_interp;
  float *disp;
  double dmean;
  Pdata_pos *moved;

  // ->> reconstruction type <<- //
  if(strcmp(rec_type, "za_displaced")==0) {
    }
  else if(strcmp(rec_type, "za_displaced_shifted")==0){
    }
  else abort();

  /* ->> get displacement field first <<- */
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  za_displacement(s, d, disp); 


  // ->> displace particles <<- //
  moved=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

  // ->> moving particles <<- //
  disp_interp=FALSE;
  move_particle(s, p, moved, disp, disp_interp);

  // ->> density <<- //
  dmean=cic_density(moved, drec, s->boxsize, s->particle_mass, s->npart, s->ngrid_xyz); 
  printf("reconstructed mean density = %lg\n", dmean);


  fftwf_free(disp);
  return;
  }
