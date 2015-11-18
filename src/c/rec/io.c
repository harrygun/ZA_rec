#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io.h"





void load_cita_simulation(char *fname, Pdata *p, int npart) {
  /* ->> load simulation data <<- */
  FILE *fd;
  int n, npt;
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
  if(npart!=npt) {printf("npart error.\n"); exit(0);}

  fread(dummy, sizeof(float), 11, fd);

  //->> loading position & velocity <<- //
  for(n=0; n<npart; n++) {
    fread(&p[n].pos[0], sizeof(float), 3, fd);
    fread(&p[n].vel[0], sizeof(float), 3, fd);
    }

  return;
  }
