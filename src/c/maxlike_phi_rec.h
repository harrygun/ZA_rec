#ifndef _H_MAX_LIKE_PHI_
#define _H_MAX_LIKE_PHI_


int phi_mlik_init(Interpar *mlik, char *fname);
double mlik_interp(Interpar *mlik, double p);



void load_stat_disp(float *disp, float *disp_phi, float *disp_model, float *phi, 
                    float *phi_model, char *fname, long long npart);


void phi_maximum_fitting(SimInfo *s, Interpar *mlik, float *phi_model, 
                        float *phi_nl, float *phi_cb, long long ngrid);

void phi_mlik_displacement(SimInfo *s, Pdata_pos *p, Interpar *mlik, float *disp, 
          float *disp_phi, float *disp_model, float *phi, float *phi_model, 
          long long ngrid, double boxsize, char *stat_disp_fname, 
          char *mlik_out_fname, int import_disp) ;


void output_maxlikelihood_data(SimInfo *s, char *fname, float *disp, float *disp_rec, 
                   float *d_rec, float *d_model, float *d_phi, long long npart) ;

#endif
