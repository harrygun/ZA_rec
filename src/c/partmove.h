#ifndef _H_PART_MOVE_
#define _H_PART_MOVE_


void move_particle(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, float *si, int s_intp);

void move_grid(SimInfo *s, Pdata_pos *moved, float *si, int s_intp);



/* general particle mover */
void move_grid_general(SimInfo *s, Pdata_pos *moved, float *si);

void general_particle_mover(Pdata_pos *p, Pdata_pos *moved, float *si, double boxsize, 
                            long long ngrid, , int s_intp);
#endif
