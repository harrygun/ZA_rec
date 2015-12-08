  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <omp.h>
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
  #include "za_reconstruction.h"
  #include "misc.h"



void za_displacement(SimInfo *s, float *d, float *disp) {
  // ->> obtain ZA displacement field from density field <<- // 
  int fft_return_type;
  float *phi, *phi_ij;

  fft_return_type=_RETURN_GRADIENT_;

  // ->> solve FFTW  <<- //
  printf("\n->> Solve Poisson equation with FFT.\n");
  poisson_solver_float(d, phi, disp, phi_ij, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R, fft_return_type);
  printf("->> FFT is Done <<- \n\n");

  return;
  }



void za_displacement_pert(SimInfo *s, float *d, float *disp) {
  // ->> obtain "perturbative-ZA" displacement field from density field <<- // 
  int fft_return_type, i, j;
  long long ip;
  float *phi, *phi_i, *phi_ij, p_ij, mat[3][3], imat[3][3], p_i[3], pc_i[3];

  fft_return_type=_RETURN_GRADIENT_HESSIAN_;

  phi_i=(float *)fftwf_malloc(sizeof(float)*s->ngrid_xyz[0]*s->ngrid_xyz[1]*s->ngrid_xyz[2]*3);
  phi_ij=(float *)fftwf_malloc(sizeof(float)*s->ngrid_xyz[0]*s->ngrid_xyz[1]*s->ngrid_xyz[2]*9);

  // ->> solve FFTW  <<- //
  printf("\n->> Solve Poisson equation with FFT.\n");
  poisson_solver_float(d, phi, phi_i, phi_ij, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R, fft_return_type);
  printf("->> FFT is Done <<- \n\n");

  //->> inverse <<- //
  #pragma omp parallel for private(ip,i,j,p_ij,mat,imat,p_i,pc_i)
  for(ip=0; ip<s->npart; ip++) {

    //-> 
    for(i=0; i<3; i++) {
      p_i[i]=ArrayAccess2D_n2(phi_i, 3, s->npart, i, ip);

      for(j=0; j<3; j++) {
	p_ij=ArrayAccess3D_n3(phi_ij, 3, 3, s->npart, i, j, ip);

        if(i==j) mat[i][j]=1.-p_ij;
	else mat[i][j]=p_ij;
        }
      }
    // ->> inverse <<- //
    mat_inv_3d_float(mat, imat);
    // ->> multiply <<- //
    mat_multiply_3d_float(imat, p_i, pc_i);

    for(i=0; i<3; i++)
      ArrayAccess2D_n2(disp, 3, s->npart, i, ip)=pc_i[i];
    }


  fftwf_free(phi_i); fftwf_free(phi_ij);
  return;
  }




void za_reconstruction(RectCtrl *rc, SimInfo *s, Pdata_pos *p, float *d, 
                         float *drec, float *d_disp, float *d_shift)   {
  // ->> Performing ZA reconstruction  <<- //
  int do_disp, do_shift, do_disp_shift, i;
  float *disp;
  double dm_disp, dm_shift;
  Pdata_pos *p_disp, *p_shift;

  // ->> reconstruction type <<- //
  if(strcmp(rc->rec_type, "displaced")==0) 
    {do_disp=TRUE;  do_shift=FALSE;  do_disp_shift=FALSE;}
  else if(strcmp(rc->rec_type, "displaced_shifted")==0)
    {do_disp=TRUE;  do_shift=TRUE;  do_disp_shift=TRUE;}
  else {abort();}

  /* ->> get displacement field first <<- */
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  if(rc->do_disp_perturb==TRUE)
    za_displacement_pert(s, d, disp); 
  else
    za_displacement(s, d, disp); 

  // ->> displace particles <<- //
  if(do_disp==TRUE) {
    p_disp=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> moving particles <<- //
    move_particle(s, p, p_disp, disp, rc->displacement_intp);

    // ->> density field from displaced particles <<- //
    dm_disp=cic_density(p_disp,d_disp,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz); 
    printf("displaced mean density delta = %lg\n", dm_disp);

    free(p_disp);
    }

  // ->> shift particles <<- //
  if(do_shift==TRUE) {
    p_shift=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));

    // ->> shifting uniform grid <<- //
    move_grid(s, p_shift, disp, rc->displacement_intp);

    // ->> density field from shifted particles <<- //
    dm_shift=cic_density(p_shift,d_shift,s->boxsize,s->particle_mass,s->npart,s->ngrid_xyz); 
    printf("shifted mean density delta = %lg\n", dm_shift);

    free(p_shift);
    }


  // ->> reconstructed density <<- //
  #pragma omp parallel for private(i)
  for(i=0; i<s->npart; i++)
    drec[i]=d_disp[i]-d_shift[i];


  fftwf_free(disp);
  return;
  }
