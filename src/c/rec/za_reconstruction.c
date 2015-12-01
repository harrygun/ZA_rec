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





void za_reconstruction(RectCtrl *rc, SimInfo *s, Pdata_pos *p, float *d, 
                         float *drec, float *d_disp, float *d_shift)   {
  // ->> Performing ZA reconstruction  <<- //
  int do_disp, do_shift, do_disp_shift;
  float *disp;
  double dm_disp, dm_shift;
  Pdata_pos *p_disp, *p_shift;

  // ->> reconstruction type <<- //
  if(strcmp(rc->rec_type, "za_displaced")==0) 
    {do_disp=TRUE;  do_shift=FALSE;  do_disp_shift=FALSE;}
  else if(strcmp(rc->rec_type, "za_displaced_shifted")==0)
    {do_disp=TRUE;  do_shift=TRUE;  do_disp_shift=TRUE;}
  else {abort();}

  /* ->> get displacement field first <<- */
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  za_displacement(s, d, disp); 

  // ->> displace particles <<- //
  if(do_disp==TRUE) {
    p_disp=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> moving particles <<- //
    move_particle(s, p, p_disp, disp, rc->displacement_intp);

    // ->> recover density <<- //
    dm_disp=cic_density(p_disp,d_disp,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz); 
    printf("displaced mean density = %lg\n", dm_disp);
    }

  // ->> shift particles <<- //
  if(do_shift==TRUE) {
    p_shift=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> moving particles <<- //
    move_particle(s, p, p_disp, disp, rc->displacement_intp);

    // ->> recover density <<- //
    dm_shift=cic_density(p_shift,d_shift,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz); 
    printf("shifted mean density = %lg\n", dm_shift);
    }


  // ->> reconstructed density <<- //


  fftwf_free(disp);
  return;
  }
