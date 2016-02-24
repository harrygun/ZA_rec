#ifndef _H_STAT_MODEL_
#define _H_STAT_MODEL_




// ->> displacement field manipulation <<- //
void load_displacement(Cospar *cp, SimInfo *s, Pdata_pos *p, float *disp, 
               float *disp_lpt, char *fname_part_init, char *fname_pid_init);

void disp_field_tranfunc_precal(SimInfo *s, Pdata_pos *p, float *d, 
                                 char *fname_part_init, char *fname_out);

//void disp_stat_separation(Cospar *cp, SimInfo *s, float *disp, float *disp_lpt, float *disp_mc, Interpar *tf);
void disp_stat_separation(Cospar *cp, float *disp, float *disp_lpt, float *disp_mc, 
                      Interpar *tf, float boxsize, long long npart, long long ngrid);


// ->> transfer function <<- //
double tk_interp(Interpar *tf, double k);


//Interpar *transfer_func_init(char *fname);
int transfer_func_init(Interpar *tf, char *fname);

void transfer_func_finalize(Interpar *tf);

//void get_stat_disp_model(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, char *stat_disp_model_type);

// -> 
int phi_mlik_init(Interpar *mlik, char *fname);
double mlik_interp(Interpar *mlik, double p);


// ->> output <<- //
void output_real_disp_field(float *disp, float *disp_lpt, size_t ngrid, char *fname_out);

void output_stat_disp_model(float *disp, float *disp_lpt, float *disp_mc, 
                      size_t ngrid, size_t ngrid_trimmed, char *fname_out);

void output_stat_disp_potential_model(float *disp, float *disp_lpt, float *disp_mc, 
                          float *div, float *phi, float *disp_phi, float *div_lpt, float *phi_lpt, 
			  float *disp_phi_lpt,  size_t ngrid, size_t ngrid_trim, char *fname_out);

#endif
