#ifdef _H_FOURIER_
#define _H_FOURIER_

//void smooth_field(float *d, double boxsize, int ngrid, int smooth_type, double smooth_R);
void smooth_field(float *d, double boxsize, int ngrid, int smooth_type, double smooth_R, Interpar *sw);

void potential_curlfree_vec(float *disp, float *div, float *phi, float *disp_phi, double boxsize, int ngrid);


void fft_gradient(float *phi, float *disp_phi, double boxsize, int ngrid);

#endif
