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



void pidtogrid(long long pid, long long ngrid, long long idx[3]) {
  // ->> return grid index of  <<- //
  //
  idx[0]=(long long)(pid/(ngrid*ngrid));
  idx[1]=(long long)(pid/ngrid)-idx[0]*ngrid;
  idx[2]=pid-ngrid*idx[1]-ngrid*ngrid*idx[0];  

  //idx[0]=((double)pid/pow((double)ngrid,2.));
  //idx[1]=((double)pid/ngrid)-idx[0]*ngrid;
  //idx[2]=pid-ngrid*idx[1]-ngrid*ngrid*idx[0];  

  return;
  }


void get_real_displacement(SimInfo *s, Pdata_pos *p, Pdata_pos *pinit, float *disp, 
                           double *drift, char *disp_calmethod, double fscale) {
  // ->> get real displacement field from simulation directly <<- //
  // ->> fscale == factor of rescaling displacement field <<- //

  long long ip, i, j, k, m, idx[3], pid, lenfac;
  float grid[3], xmin[3], xmax[3], maxdisp[3], dx;
  double disp_, fac;

  lenfac=12;
  fac=0.95;

  // ->> box boundary <<- //
  dx=s->boxsize/(float)s->ngrid;

  printf("max disp:  ");
  for(i=0; i<3; i++){
    xmin[i]=0.5;//drift[2-i];
    xmax[i]=xmin[i]+(s->ngrid)*dx;

    //maxdisp[i]=min(s->boxsize-drift[i]*dx, s->boxsize+drift[i]*dx);
    maxdisp[i]=min(s->boxsize-drift[2-i]*dx, s->boxsize+drift[2-i]*dx);
    printf("%lg  (%lg)  ", maxdisp[i], maxdisp[i]*fac);
    }
  printf("\n"); fflush(stdout);


  double *pmin, *pmax, *dpart;
  pmin=(double *)malloc(3*sizeof(double));
  pmax=(double *)malloc(3*sizeof(double));
  dpart=(double *)malloc(3*sizeof(double));

  get_particle_boundary(p, s->boxsize, s->npart, s->ngrid_xyz, pmin, pmax, dpart);
  printf("->> xmin %f ymin %f zmin %f\n",pmin[0], pmin[1], pmin[2]);
  printf("->> xmax %f ymax %f zmax %f\n",pmax[0], pmax[1], pmax[2]);


  // ->>
  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
    printf("direct subtraction for real_displacement.\n"); fflush(stdout);

    #ifdef _OMP_
    #pragma omp parallel for private(ip,i)
    #endif
    for(ip=0; ip<s->npart; ip++) {
      for(i=0; i<3; i++) {
        ArrayAccess2D_n2(disp, 3, s->npart, i, ip)=(p[ip].pos[i]-pinit[ip].pos[i])*fscale;
        }
      }
    }
  else if( strcmp(disp_calmethod, "grid_wise")==0 ){
    printf("grid-wise calculation for real_displacement.\n"); fflush(stdout);

    // ->> re-arrange data into grid <<- //
    #ifdef _OMP_
    #pragma omp parallel for private(i,j,k,m,ip,grid,disp_)
    #endif
    for(i=0; i<s->ngrid; i++)
      for(j=0; j<s->ngrid; j++)
        for(k=0; k<s->ngrid; k++) {

          // ->> grid index <<- //
          ip=MemIdx3D(s->ngrid, i, j, k);

          grid[0]=xmin[0]+i*dx;
          grid[1]=xmin[1]+j*dx;
          grid[2]=xmin[2]+k*dx;

          for(m=0; m<3; m++){
            //disp_=p[ip].pos[2-m]-grid[m];
            ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=(p[ip].pos[2-m]-grid[m])*fscale;
            //ArrayAccess2D_n2(disp, 3, s->npart, m, ip)=(p[ip].pos[2-m]-drift[2-m]-grid[m])*fscale;
            }

          }
    }
  else if( strcmp(disp_calmethod, "grid_PID")==0 ){
    printf("grid-PID calculation for real_displacement.\n"); fflush(stdout);

    #ifdef _OMP_
    #pragma omp parallel for private(ip,i,idx,grid,pid,disp_)
    #endif
    for(ip=0; ip<s->npart; ip++) {
      //->> get the index <<- //
      pid=p[ip].pid-1;
      pidtogrid(pid, s->ngrid, idx);

      for(i=0; i<3; i++) {
        grid[i]=xmin[i]+idx[i]*dx;
        //disp_=(p[ip].pos[2-i]-pmin[2-i])-grid[i];
        disp_=(p[ip].pos[2-i]-drift[2-i])-grid[i];


        if(disp_<=-maxdisp[i]*fac)
          disp_+=maxdisp[i]*fac;
        if(disp_>=maxdisp[i]*fac) 
          disp_-=maxdisp[i]*fac;

        ArrayAccess2D_n2(disp, 3, s->npart, i, pid)=disp_*fscale;
        }
      }
    }
  else if( strcmp(disp_calmethod, "init_pos_PID")==0 ){
    printf("init_pos-PID calculation for real_displacement.\n"); fflush(stdout);

    #ifdef _OMP_
    #pragma omp parallel for private(ip,i,j,idx,pid,disp_)
    #endif
    for(ip=0; ip<s->npart; ip++) {
      //->> get the index <<- //
      pid=p[ip].pid-1;
      pidtogrid(pid, s->ngrid, idx);

      for(i=0; i<3; i++) {
        //grid[i]=xmin[i]+idx[i]*dx;
        //disp_=(p[ip].pos[2-i]-pmin[2-i])-pinit[ip].pos[2-i];

        disp_=p[ip].pos[2-i]-pinit[pid].pos[2-i];

        /*
        for(j=0; j<lenfac; j++){
          fac=1-j*0.005;

          if(disp_<=-maxdisp[i]*fac)
            disp_+=maxdisp[i]*fac;
          if(disp_>=maxdisp[i]*fac) 
            disp_-=maxdisp[i]*fac;

	  }
	*/
        fac=0.95;
        if(fabs(disp_)>=maxdisp[i]*fac) disp_=0.;
	

        ArrayAccess2D_n2(disp, 3, s->npart, i, pid)=disp_*fscale;
        }
      }


    }

  else {abort();}

  return;
  }





void get_model_displacement(SimInfo *s, Pdata_pos *p, float *d, float *disp, char *fname_part_init, char *modeldisp_type){

  int model_from_init_pos=TRUE;

  Pdata_pos *pinit;
  double dmean;

  if (model_from_init_pos==TRUE){
    pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
    load_cita_simulation_position(fname_part_init, pinit, s->npart);
      
    // ->> estimate the density field from initial field <<- //
    dmean=cic_density(pinit, d, s->boxsize, s->particle_mass, s->npart, s->ngrid_xyz, s); 
    }


  if(strcmp(modeldisp_type, "ZA")==0 ) {
    za_displacement(s, d, disp);
    }
  else if(strcmp(modeldisp_type, "2LPT")==0 ) {
    displacement_2lpt(s, d, disp);
    }
  else abort();


  if (model_from_init_pos==TRUE){
    free(pinit); }

  return;
  }



/* ->> some manupulations of displacement fields <<- */
/*
void disp_potential(SimInfo *s, float *disp, float *disp_phi){
  int i; 
  float *phi, *phi_i, *phi_ij;


  for(i=0; i<3; i++) {

    //ArrayAccess2D_n2(disp, 3, s->npart, i, 0);

    poisson_solver_float(ArrayAccess2D_n2(disp, 3, s->npart, i, 0), phi, disp, phi_ij, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R, fft_return_type);
    }


  return;
  }

*/



