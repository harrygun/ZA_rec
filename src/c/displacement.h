#ifdef _H_DISPLACEMENT_
#define _H_DISPLACEMENT_


void get_real_displacement(SimInfo *s, Pdata_pos *p, float *disp, 
                           char *fname_part_init, char *disp_calmethod);

//void get_model_displacement(SimInfo *s, Pdata_pos *p, char *model_disp_type);
void get_model_displacement(SimInfo *s, Pdata_pos *p, float *d, float *disp, char *fname_part_init, char *model_disp_type);


#endif
