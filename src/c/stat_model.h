#ifndef _H_STAT_MODEL_
#define _H_STAT_MODEL_

Interpar *transfer_func_init(char *fname);



// ->> displacement field manipulation <<- //
void load_displacement(Cospar *cp, SimInfo *s, Pdata_pos *p, float *disp, 
                         float *disp_lpt, char *fname_part_init);

void disp_field_tranfunc_precal(SimInfo *s, Pdata_pos *p, float *d, 
                                 char *fname_part_init, char *fname_out);

void disp_stat_separation(Cospar *cp, SimInfo *s, float *disp, float *disp_lpt, float *disp_mc, Interpar *tf);
void output_stat_disp_model(float *disp, float *disp_lpt, float *disp_mc, char *fname_out);


// ->> transfer function <<- //
double tk_interp(Interpar *tf, double k);
Interpar *transfer_func_init(char *fname);
void transfer_func_finalize(Interpar *tf);

//void get_stat_disp_model(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, char *stat_disp_model_type);

#endif
