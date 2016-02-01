#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io.h"





void load_cita_simulation(char *fname, Pdata *p, long long npart) {
  /* ->> load simulation data <<- */
  FILE *fd;
  long long n, npt;
  float dummy[11];
    
  // ->> start to read <<- //
  if(!(fd=fopen(fname, "r"))) {
    printf("can't open file `%s`\n", fname);
    exit(0);
    }
 
  printf("reading `%s' ...\n", fname);
  fflush(stdout);

  fread(&npt, 4, 1, fd);
  printf("npt=%d, npart=%d\n", npt, npart);
  if((int)npart!=(int)npt) {printf("npart error.\n"); exit(0);}

  fread(dummy, sizeof(float), 11, fd);

  //->> loading position & velocity <<- //
  for(n=0; n<npart; n++) {
    fread(&p[n].pos[0], sizeof(float), 3, fd);
    fread(&p[n].vel[0], sizeof(float), 3, fd);
    }

  return;
  }



void load_cita_simulation_position(char *fname, Pdata_pos *p, long long npart) {
  /* ->> load only the position of simulation data <<- */
  FILE *fd;
  long long n, npt;
  float dummy[11], vel[3];
    
  // ->> start to read <<- //
  if(!(fd=fopen(fname, "r"))) {
    printf("can't open file `%s`\n", fname);
    exit(0);
    }
 
  printf("reading `%s' ...\n", fname);
  fflush(stdout);

  fread(&npt, 4, 1, fd);
  printf("npt=%d, npart=%d\n", npt, npart);
  if((int)npart!=(int)npt) {printf("npart error.\n"); exit(0);}

  fread(dummy, sizeof(float), 11, fd);

  //->> loading position & velocity <<- //
  for(n=0; n<npart; n++) {
    fread(&p[n].pos[0], sizeof(float), 3, fd);
    fread(vel, sizeof(float), 3, fd);
    }

  return;
  }



void load_scalar_map(char *fname, float *m, long long ngrid, char *dtype){
  int size;
  FILE *fp=fopen(fname, "r");

  if(strcmp(dtype, "float")==0 ) 
    size=sizeof(float);
  else if(strcmp(dtype, "double")==0 ) 
    size=sizeof(double);
  else
    abort();

  fread(m, size, ngrid*ngrid*ngrid, fp);
  fclose(fp);
  return;
  }
