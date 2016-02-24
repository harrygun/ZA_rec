#ifndef _H_IO_
#define _H_IO_

typedef struct particle_data {
  float pos[3];
  float vel[3];
  //float Mass;
  //int Type;
  //float Rho, U, Temp, Ne;
  long long pid;
  } Pdata;


typedef struct particle_pos_data {
  float pos[3];
  long long pid;
  } Pdata_pos;

void load_simulation_offset(char *fname, double *offset_f, double *offset_i);


void load_cita_simulation(char *fname, Pdata *p, long long npart);
void load_cita_simulation_position(char *fname, Pdata_pos *p, long long npart);

void load_cita_simulation_pid(char *fname, Pdata_pos *p, long long npart);

void load_cita_simulation_position_pid(char *fname_pos, char *fname_pid, Pdata_pos *p, long long npart);


void load_scalar_map(char *fname, float *m, long long ngrid, char *dtype);
#endif
