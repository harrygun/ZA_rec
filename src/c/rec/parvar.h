#ifndef _H_PARVAR_
#define _H_PARVAR_

  #define _TOPHAT_SMOOTH_  20121201
  #define _GAUSSIAN_SMOOTH_ 20150228

  // ->> for equal-length cubic array <<- //
  #define ArrayAccess3D(a, n, i, j, k) ((a)[(i)*(n)*(n)+(j)*(n)+(k)])

  // ->> non-equal-length cubic array <<- //
  #define ArrayAccess3D_n3(a, n1, n2, n3, i, j, k) ((a)[ ((n2)*(i)+(j))*(n3)+(k) ])


  #define ArrayAccess5D_n5(a, n1, n2, n3, n4, n5, i1, i2, i3, i4, i5) ((a)[ i5+n5*(i4+n4*(i3+n3*(i2+n2*i1))) ])
  #define ArrayAccess4D_n4(a, n1, n2, n3, n4, i1, i2, i3, i4) ((a)[ i4+n4*(i3+n3*(i2+n2*i1)) ])
   
  #define ArrayAccess2D_n2(a, n1, n2, i1, i2) ((a)[ i2+n2*i1 ])
   

#endif
