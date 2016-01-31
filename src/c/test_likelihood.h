#ifndef _H_TEST_LIKELIHOOD_
#define _H_TEST_LIKELIHOOD_


void test_displacement(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                       char *fname_out);



void test_disp_direct_cal(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *fname_out);
void test_disp_vel_comp(SimInfo *s, Pdata_pos *p, float *d, char *fname_part_init, 
                          char *fname_out);
#endif
