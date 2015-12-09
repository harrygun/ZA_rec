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
  #include "poisson.h"
  #include "misc.h"
  #include "backward_displacement.h"

  #ifdef _OMP_
  #include <omp.h>
  #endif





void reconstruction_partmover(RectCtrl *rc, SimInfo *s, Pdata_pos *p, float *d, 
                                 float *drec, float *d_disp, float *d_shift)   {
  // ->> Performing ZA reconstruction  <<- //
  int i, do_disp, do_shift, do_disp_shift, lpt_order;
  float *disp;
  double dm_disp, dm_shift;
  Pdata_pos *p_disp, *p_shift;

  // ->> reconstruction type <<- //
  if(strcmp(rc->rec_type, "displaced")==0) 
    {do_disp=TRUE;  do_shift=FALSE;  do_disp_shift=FALSE;}
  else if(strcmp(rc->rec_type, "displaced_shifted")==0)
    {do_disp=TRUE;  do_shift=TRUE;  do_disp_shift=TRUE;}
  else {abort();}

  // ->> forward vs. backward displacement <<- //
  if(strcmp(rc->displacement_type, "backward_displacement")!=0) {
    printf("only support backward modeling now.\n");
    fflush(stdout); abort();
    }

  // ->> LPT order <<- //
  if(strcmp(rc->displacement_order, "1LPT")==0) 
    {lpt_order=1;}
  else if(strcmp(rc->displacement_order, "2LPT")==0) 
    {lpt_order=2;}
  else {lpt_order=-1; abort();}


  /* ->> obtain the displacement field from density field <<- */
  // ->> memory allocation first <<- //
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  if(lpt_order==1){
    // ->> Zel'dovich Approximation <<- //

    if(rc->do_disp_perturb==TRUE) {
      printf("perform perturbed-ZA reconstruction.\n"); fflush(stdout);
      za_displacement_pert(s, d, disp); 
      }
    else{
      printf("perform unperturbed-ZA reconstruction.\n"); fflush(stdout);
      za_displacement(s, d, disp); 
      }
    }
  else if(lpt_order==2) {
    // ->> 2-LPT <<- //
    if(rc->do_disp_perturb==TRUE) {
      printf("`perturbative' 2LPT NOT supported yet.");
      fflush(stdout); abort();
      }
    else
      displacement_2lpt(s, d, disp); 
    }


  /* ->> moving particles <<- */
  if(do_disp==TRUE) {
    p_disp=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> moving particles <<- //
    move_particle(s, p, p_disp, disp, rc->displacement_intp);

    // ->> density field from displaced particles <<- //
    dm_disp=cic_density(p_disp,d_disp,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz, NULL); 
    printf("displaced mean density delta = %lg\n", dm_disp);

    free(p_disp);
    }

  // ->> shift particles <<- //
  if(do_shift==TRUE) {
    p_shift=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> shifting uniform grid <<- //
    move_grid(s, p_shift, disp, rc->displacement_intp);

    // ->> density field from shifted particles <<- //
    dm_shift=cic_density(p_shift,d_shift,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz, NULL); 
    printf("shifted mean density delta = %lg\n", dm_shift);

    free(p_shift);
    }


  // ->> reconstructed density <<- //
  #ifdef _OMP_
  #pragma omp parallel for private(i)
  #endif
  for(i=0; i<s->npart; i++)
    drec[i]=d_disp[i]-d_shift[i];


  fftwf_free(disp);
  return;
  }
