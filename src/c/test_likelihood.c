  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <iniparser.h>

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
  #include "reconstruction_partmoving.h"
  #include "fourier.h"
  #include "backward_displacement.h"

  #include "stat_model.h"



#ifdef _OMP_
  #include <omp.h>
#endif




void test_displacement(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                       char *fname_out) {
  // ->> try to build the statistical model from real displacement <<- //
  //get_stat_disp_model(s, p, d, fname_part_init, NULL);
  int i;
  float *disp, *disp_model;
  Pdata_pos *pinit;

  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  disp_model=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);

  // ->> obtain real displacement <<- //
  char *disp_calmethod="grid_wise";

  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
    // ->> load initial position <<- //
    pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata));
    load_cita_simulation_position(fname_part_init, pinit, s->npart);
    }

  get_real_displacement(s, p, pinit, disp, disp_calmethod, 1.);


  // ->> smooth the field <<- //
  for(i=0; i<3; i++) {
    smooth_field(&disp[i*s->ngrid*s->ngrid*s->ngrid], s->boxsize, 
                       s->ngrid, s->smooth_type_flag, s->smooth_R);
    }


  // ->> obtain model displacement <<- //
  char *model_type="ZA"; //"2LPT";  //"ZA";
  get_model_displacement(s, p, d, disp_model, fname_part_init, model_type);

  // ->> construct model <<- //
  FILE *fp=fopen(fname_out, "wb");

  fwrite(disp, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(disp_model, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(d, sizeof(float), s->ngrid*s->ngrid*s->ngrid, fp);

  fclose(fp);


  // ->> free <<- //
  free(disp); free(disp_model);
  if(strcmp(disp_calmethod, "direct_subtraction")==0 ) {
    free(pinit); }

  return;
  }




void test_disp_direct_cal(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *fname_out) {
  // ->> direct estimation of displacement field <<- //
  int i;
  float *disp_init, *disp;
  double dmean;
  disp_init=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  disp=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);


  // ->> get displacement field differently <<- //
  Pdata_pos *pinit=(Pdata_pos *)malloc(s->npart*sizeof(Pdata_pos));
  load_cita_simulation_position(fname_part_init, pinit, s->npart);

  // ->>  get initial displacement <<- //
  char *disp_calmethod="grid_wise";
  get_real_displacement(s, pinit, pinit, disp_init, disp_calmethod, 1.);

  // ->>  get final displacement <<- //
  get_real_displacement(s, p, pinit, disp, disp_calmethod, 1.);
   

  // ->> obtain model displacement <<- //
  dmean=cic_density(p, d, s->boxsize, s->particle_mass, s->npart, s->ngrid_xyz, s); 
  //za_displacement(s, d, disp);

  // ->> construct model <<- //
  FILE *fp=fopen(fname_out, "wb");

  fwrite(disp_init, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(disp, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(d, sizeof(float), s->ngrid*s->ngrid*s->ngrid, fp);

  fclose(fp);

  // ->> free <<- //
  free(disp_init); free(disp);
  free(pinit);

  return;
  }





// ->>  <<- //
void test_disp_vel_comp(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                        char *fname_out) {
  int i;
  long long ip;
  float *disp_init, *vel, *disp_model;
  double dmean;
  disp_init=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  vel=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);


  // ->> get displacement field differently <<- //
  Pdata *pinit=(Pdata *)malloc(s->npart*sizeof(Pdata));
  load_cita_simulation(fname_part_init, pinit, s->npart);

  // ->>
  #ifdef _OMP_
  #pragma omp parallel for private(ip,i)
  #endif
  for(ip=0; ip<s->npart; ip++) {
    for(i=0; i<3; i++) {
      ArrayAccess2D_n2(vel, 3, s->npart, i, ip)=pinit[ip].vel[i];
      }
    }

  // ->> obtain model displacement <<- //
  Pdata_pos *pos_init=(Pdata_pos *)malloc(s->npart*sizeof(Pdata_pos));
  load_cita_simulation_position(fname_part_init, pos_init, s->npart);
  dmean=cic_density(pos_init, d, s->boxsize, s->particle_mass, s->npart, s->ngrid_xyz, s); 


  // ->>  grid-wise displacement field <<- //
  char *disp_calmethod="grid_wise";
  get_real_displacement(s, p, pos_init, disp_init, disp_calmethod, 1.);


  // ->> ZA <<- //
  //disp_model=(float *)fftwf_malloc(sizeof(float)*s->ngrid*s->ngrid*s->ngrid*3);
  //s->smooth_R=0.;
  //za_displacement(s, d, disp_model);

  // ->> construct model <<- //
  FILE *fp=fopen(fname_out, "wb");

  fwrite(disp_init, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  //fwrite(disp_model, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(vel, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  fwrite(d, sizeof(float), s->ngrid*s->ngrid*s->ngrid, fp);

  fclose(fp);

  // ->> free <<- //
  free(disp_init); free(vel);
  free(pinit);  free(pos_init);

  return;
  }



