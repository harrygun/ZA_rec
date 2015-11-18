#ifndef _H_IO_
#define _H_IO_

typedef struct cita_particle_data
{
  float pos[3];
  float vel[3];
  //float Mass;
  //int Type;
  //float Rho, U, Temp, Ne;
} Pdata;


void load_cita_simulation(char *fname, Pdata *p, int npart);

#endif
