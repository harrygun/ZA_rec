#ifndef _H_CIC_
#define _H_CIC_


double get_rhom(Cospar *cp, double z);
double part_mass(Cospar *cp, double z, double boxsize, int ngrid);


//double cic_density(Pdata_pos *p, double ***d, double boxsize, double mass, int npart, int ngrid[3]);
//double cic_density(Pdata_pos *p, double *d, double boxsize, double mass, int npart, int ngrid[3]);
double cic_density(Pdata_pos *p, float *d, double boxsize, double mass, int npart, int ngrid[3]);

#endif
