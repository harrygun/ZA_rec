#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "readfile.h"
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


  fclose(fd);
  return;
  }

void load_cita_simulation_pid(char *fname, Pdata_pos *p, long long npart){
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
    fread(&p[n].pid, sizeof(long long), 1, fd); 
    }

  fclose(fd);
  return;
  }



// ->> load position & PID together <<- //
void load_cita_simulation_position_pid(char *fname_pos, char *fname_pid, Pdata_pos *p, long long npart){

  // ->> load position first <<- //
  load_cita_simulation_position(fname_pos, p, npart);

  // ->> then PID <<- //
  load_cita_simulation_pid(fname_pid, p, npart);

  return;
  }


void load_simulation_offset(char *fname, double *offset){
  int i;
  double z, off_i[3], off_f[3];
  FILE *fp=fopen(fname, "r");  

  fscanf(fp, "%lg ", &z);
  printf("read simulation offset file %s at redshift %lg.\n", fname, z);
  fflush(stdout);

  for(i=0; i<3; i++) {
    fscanf(fp, "%lg ", &offset[i]);
    printf("%lg ", offset[i]);
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
