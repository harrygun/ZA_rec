#ifndef _H_PARVAR_
#define _H_PARVAR_

  #define _TOPHAT_SMOOTH_  20121201
  #define _GAUSSIAN_SMOOTH_ 20150228


  typedef float itreal;
  //typedef double itreal;



typedef struct simulation_info{
  // ->> simulation information <<- //
  double boxsize, smooth_R, particle_mass, xmin[3], xmax[3];
  int smooth_type_flag;
  long long ngrid, ngrid_xyz[3], npart;
  char *smooth_type, *test_fname;
  } SimInfo;

typedef struct rect_controller{
  //->> reconstruction controller <<- 
  int do_rect, displacement_intp;
  char *rec_type, rec_fname[200];   //, *rec_fname;
  }RectCtrl;


  // ->>
  #define MemIdx3D(n, i1, i2, i3)  (i3+n*(i2+n*i1))
  #define MemIdx3D_n3(n1, n2, n3, i1, i2, i3)  (i3+n3*(i2+n2*i1))


  // ->> for equal-length cubic array <<- //
  #define ArrayAccess3D(a, n, i, j, k) ((a)[(i)*(n)*(n)+(j)*(n)+(k)])

  // ->> non-equal-length cubic array <<- //
  #define ArrayAccess3D_n3(a, n1, n2, n3, i, j, k) ((a)[ ((n2)*(i)+(j))*(n3)+(k) ])

  #define ArrayAccess5D_n5(a, n1, n2, n3, n4, n5, i1, i2, i3, i4, i5) ((a)[ i5+n5*(i4+n4*(i3+n3*(i2+n2*i1))) ])

  #define ArrayAccess4D_n4(a, n1, n2, n3, n4, i1, i2, i3, i4) ((a)[ i4+n4*(i3+n3*(i2+n2*i1)) ])
   
  #define ArrayAccess2D_n2(a, n1, n2, i1, i2) (a)[ i2+n2*i1 ]

  #define ArrayAccess2D_n2_list(a, n1, n2, nlist, i1, i2) (a)[ (i2+n2*i1)*nlist ]
   


#endif
