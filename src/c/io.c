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
void load_cita_simulation_position_pid(char *fname_pos, char *fname_pid, 
                                       Pdata_pos *p, long long npart){

  // ->> load position first <<- //
  load_cita_simulation_position(fname_pos, p, npart);

  // ->> then PID <<- //
  load_cita_simulation_pid(fname_pid, p, npart);

  return;
  }


void load_simulation_offset(char *fname, double *offset_f, double *offset_i){
  int i;
  double zi, zf, off_i[3], off_f[3];
  FILE *fp=fopen(fname, "r");  

  fscanf(fp, "%lg ", &zf);
  for(i=0; i<3; i++) {
    fscanf(fp, "%lg ", &offset_f[i]); }
  fscanf(fp, "\n");

  fscanf(fp, "%lg ", &zi);
  for(i=0; i<3; i++) {
    fscanf(fp, "%lg ", &offset_i[i]); }

  printf("read simulation offset file %s at redshift %lg (zi=%lg).\n", fname, zf, zi);
  for(i=0; i<3; i++) 
    printf("%lg   ", offset_i[i]);
  printf("\n");

  for(i=0; i<3; i++) 
    printf("%lg   ", offset_f[i]);
  printf("\n");

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




void cp_pdata_info(SimInfo *s, Pdata_pos *p) {
  /* ->> copy information of simulation <<- */

  p->npart=s->npart;   
  p->ngrid=s->ngrid; 
  p->bsize=s->boxsize;
  
  p->ng_trim=s->ng_trim; 
  p->ntrim=s->ntrim;  
  p->npart_trim=s->npart_trim; 
  p->bsize_trim=s->bsize_trim;

  if(p->ntrim>0)  {p->do_trim=TRUE;}
  else {p->do_trim=FALSE;}

  return;
  }




void pdata_access_trim(Pdata_pos *p, Pdata_pos r, size_t ip, size_t i){
  /*->> call this routine if need to trim data, store returned in 'r' <<- */
  size_t ip_old; 
  if(p->do_trim==TRUE)  {

    }
  else {

    }
  return;
  }




void pdata_trim(Pdata_pos *p, Pdata_pos *r){
  /*->> trim the pdata and stored in 'r' <<- */
  long long ip, pid, i, idx[3];

  if(p->do_trim!=TRUE) {
    r=p; return; }

  #ifdef _OMP_
  #pragma omp parallel for private(ip,i)
  #endif
  for(ip=0; ip<p->npart_trim; ip++) {

    //pid=
    //pidtogrid(pid, long long ngrid, long long idx[3]);

    for(i=0; i<3; i++) {

      }


    }

  return;
  }


