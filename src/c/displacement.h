#ifdef _H_DISPLACEMENT_
#define _H_DISPLACEMENT_



void get_real_displacement(SimInfo *s, Pdata_pos *p, Pdata_pos *pinit, float *disp, 
                           char *disp_calmethod, double fscale);

//void get_model_displacement(SimInfo *s, Pdata_pos *p, char *model_disp_type);
void get_model_displacement(SimInfo *s, Pdata_pos *p, float *d, float *disp, char *fname_part_init, char *model_disp_type);


#endif
