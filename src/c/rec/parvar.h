#ifndef _H_PARVAR_
#define _H_PARVAR_

  #define _TOPHAT_SMOOTH_  20121201
  #define _GAUSSIAN_SMOOTH_ 20150228

  // ->> for equal-length cubic array <<- //
  #define ArrayAccess3D(a, n, i, j, k) ((a)[(i)*(n)*(n)+(j)*(n)+(k)])

  // ->> non-equal-length cubic array <<- //
  #define ArrayAccess3D_n3(a, n1, n2, n3, i, j, k) ((a)[ ((n2)*(i)+(j))*(n3)+(k) ])

#endif
