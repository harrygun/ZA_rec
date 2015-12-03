#ifndef _H_IO_
#define _H_IO_

typedef struct particle_data {
  float pos[3];
  float vel[3];
  //float Mass;
  //int Type;
  //float Rho, U, Temp, Ne;
  } Pdata;


typedef struct particle_pos_data {
  float pos[3];
  } Pdata_pos;


void load_cita_simulation(char *fname, Pdata *p, int npart);
void load_cita_simulation_position(char *fname, Pdata_pos *p, int npart);

void load_scalar_map(char *fname, float *m, int ngrid, char *dtype);
#endif
