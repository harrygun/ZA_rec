  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "const.h"
  #include "parvar.h"
  #include "io.h"
  #include "cic.h"

  #include "displacement.h"

  #ifdef _OMP_
  #include <omp.h>
  #endif


float trilinear(float pos[3], float dp[3], float v[2][2][2]){

  float c00, c01, c10, c11, c0, c1, c;
  
  c00=v[0][0][0]*(1.-dp[0])+v[1][0][0]*dp[0];
  c10=v[0][1][0]*(1.-dp[0])+v[1][1][0]*dp[0];
  c01=v[0][0][1]*(1.-dp[0])+v[1][0][1]*dp[0];
  c11=v[0][1][1]*(1.-dp[0])+v[1][1][1]*dp[0];
  
  c0=c00*(1.-dp[1]) + c10*dp[1];
  c1=c01*(1.-dp[1]) + c11*dp[1];
  
  c=c0*(1.-dp[2])+c1*dp[2];
  
  return c;
  }





void move_particle(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, float *si, int s_intp){
  //-> move particles <<- //
  long long i, ip, m, n, l, m1, n1, l1, idx[3];
  float moved_pos, dsi, dp[3], v[2][2][2], xmin, xmax, dx;
  float disp__;

  // ->> box boundary <<- //
  xmin=0.; xmax=s->boxsize;
  dx=(xmax-xmin)/(float)s->ngrid;
       
  printf("->> displaceing particles, xmin=%f, xmax=%f, dx=%f\n", xmin, xmax, dx);
  //printf("->> displaceing particles ... \n");
  //

  //FILE *fp=fopen("si.dat", "wb");
  //fwrite(si, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);
  //fclose(fp);


  #ifdef _OMP_
  #pragma omp parallel for private(ip,i,m,n,l,m1,n1,l1,idx,moved_pos,dsi,dp,v)
  #endif
  for(ip=0; ip<s->npart; ip++) {

    // ->> position index <<- //
    for(i=0; i<3; i++){
      idx[i]=(long long)(((double)p[ip].pos[i]-s->pmin[i])/s->dpart[i]); 
      if(idx[i]==s->ngrid) {idx[i]=0;}

      if((idx[i]<0)||(idx[i]>s->ngrid-1)) {
        printf("move_part idx error: %d (ip=%d, i=%d), ", idx[i], ip, i); 
	printf(" %f, %f, %f, %f, foat(idx)=%f\n", p[ip].pos[i], (float)s->pmin[i], (float)s->pmax[i], 
	       (float)s->dpart[i], (p[ip].pos[i]-(float)s->pmin[i])/(float)s->dpart[i]);
	fflush(stdout);}
      }

    // ->> do not interpolate <<- //
    if(s_intp==FALSE){
      //->> particle-moving DO NOT interpolation.\n");

      for(i=0; i<3; i++){
        moved_pos=p[ip].pos[i]+ArrayAccess4D_n4(si, 3, s->ngrid, s->ngrid, s->ngrid, i, idx[0], idx[1], idx[2]);
	

        // ->> periodic boundary condition <<- //
	if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
	else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
	else{
	  moved[ip].pos[i]=moved_pos; }
	}
      }

    // ->> interpolate shift field onto particle position <<- //
    else if(s_intp==TRUE){  
      // ->> ("particle-moving interpolation.\n") <<- //
 
      // ->> distance to grid <<- //
      for(i=0; i<3; i++) {
        dp[i]=(float)((double)p[ip].pos[i]-s->pmin[i]-idx[i]*s->dpart[i]); 
	if(dp[i]>=(float)(s->pmax[i]-s->pmin[i])) dp[i]=0.; 

	if((dp[i]<0)||(dp[i]>(float)s->dpart[i]))  {
	  printf("dp error, dp=%f (pos=%f, pos_grid=%f, idx=%d, dpart=%f)\n", dp[i], p[ip].pos[i], 
	         idx[i]*(float)s->dpart[i], idx[i], (float)s->dpart[i]);  
	  fflush(stdout);
	  }

	}

      for(i=0; i<3; i++)  {

        // ->> define vertices for 3D interpolation <<- //
        for(l=0; l<2; l++){
	  l1=l+idx[0]; 
	  if(l1>=s->ngrid)  l1=l1-s->ngrid;

          for(m=0; m<2; m++) {
	    m1=m+idx[1]; 
	    if(m1>=s->ngrid)  m1=m1-s->ngrid;

            for(n=0; n<2; n++) {
	      n1=n+idx[2]; 
	      if(n1>=s->ngrid)  n1=n1-s->ngrid;


              v[l][m][n]=ArrayAccess4D_n4(si, 3, s->ngrid, s->ngrid, s->ngrid, i, l1, m1, n1);
	      }
	    }
	  }

        // ->> tri-linear interpolation <<- //
        dsi=trilinear(p[ip].pos, dp, v);
        moved_pos=p[ip].pos[i]+dsi;

	if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
	else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
	else{
	  moved[ip].pos[i]=moved_pos; }
        }
      }

    else{abort();}

    }


  //#define _MOVED_PART_OUTPUT_
  #ifdef _MOVED_PART_OUTPUT_
  FILE *fp=fopen(s->test_fname, "wb");
  //fp=fopen("si_after.dat", "wb");

  for(ip=0; ip<s->npart; ip++) 
    fwrite(moved[ip].pos, sizeof(float), 3, fp);
 
  fwrite(si, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);

  fclose(fp);
  #endif


  printf("->> particles displacement is done.\n");
  return;
  }



void move_grid(SimInfo *s, Pdata_pos *moved, float *si, int s_intp){
  //->> grid moving, no need to generate the grid <<- //
  long long i, j, k, m, ip, idx[3];
  float grid[3], pos[3], xmin, xmax, dx, moved_pos;

  //->> do not need to interpolate, on grid already <<- //
  s_intp=FALSE;
       
  // ->> box boundary <<- //
  xmin=0.; xmax=s->boxsize;
  dx=(xmax-xmin)/(float)s->ngrid;

  printf("\n->> shifting uniform grid ...\n");
  

  FILE *fp=fopen("/mnt/scratch-lustre/xwang/data/baorec/random/random_ng512_bx1024.dat", "rb");
  float *rand=(float *)malloc(sizeof(float)*s->npart*3);
  fread(rand, sizeof(float), s->npart*3, fp);
  fclose(fp);


  #ifdef _OMP_
  #pragma omp parallel for private(i,ip,pos,idx,moved_pos)
  #endif
  for(ip=0; ip<s->npart; ip++) {

    // ->> random particle 
    pos[0]= rand[ip];            
    pos[1]= rand[s->npart+ip];   
    pos[2]= rand[2*s->npart+ip]; 


    // ->> position index <<- //
    for(i=0; i<3; i++){

      idx[i]=(long long)(((double)pos[i])/s->dpart[i]); 
      if(idx[i]==s->ngrid) {idx[i]=0;}

      if((idx[i]<0)||(idx[i]>s->ngrid-1)) {
        printf("move_part idx error: %d (ip=%d, i=%d), ", idx[i], ip, i); 
	printf(" %f, %f, %f, %f, foat(idx)=%f\n", pos[i], (float)s->pmin[i], (float)s->pmax[i], 
	       (float)s->dpart[i], (pos[i]-(float)s->pmin[i])/(float)s->dpart[i]);
	fflush(stdout);}
      }

      for(i=0; i<3; i++)  {
        moved_pos=pos[i]+ArrayAccess4D_n4(si, 3, s->ngrid, s->ngrid, s->ngrid, i, idx[0], idx[1], idx[2]);
      
        // ->> periodic boundary condition <<- //
        if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
        else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
        else{
           moved[ip].pos[i]=moved_pos; }
         }

    }


  /*
  #ifdef _OMP_
  #pragma omp parallel for private(i,j,k,m,ip,grid,moved_pos)
  #endif
  for(i=0; i<s->ngrid; i++)
    for(j=0; j<s->ngrid; j++)
      for(k=0; k<s->ngrid; k++) {

        // ->> do not interpolate <<- //
        if(s_intp==FALSE){
          // ->> grid index <<- //
          ip=MemIdx3D(s->ngrid, i, j, k);

	  grid[0]= rand[ip];            //xmin+i*dx;
	  grid[1]= rand[s->npart+ip];   //xmin+j*dx;
	  grid[2]= rand[2*s->npart+ip]; //xmin+k*dx;

	  for(m=0; m<3; m++)  {
            moved_pos=grid[m]+ArrayAccess2D_n2(si, 3, s->npart, m, ip);

            // ->> periodic boundary condition <<- //
	    if(moved_pos<0){
              moved[ip].pos[m]=moved_pos+xmax; }
	    else if(moved_pos>=xmax){
              moved[ip].pos[m]=moved_pos-xmax; }
	    else{
	      moved[ip].pos[m]=moved_pos; }
            }

          }

        // ->> interpolate shift field onto particle position <<- //
        else if(s_intp==TRUE){  
          printf("grid-moving interpolation NOT supported yet.\n");
          abort();
          }
        else{abort();}
      }
  */




  printf("->> shifting grid is done.\n");
  return;
  }



/* ->> General particle moving routines <<- */
void general_particle_mover(SimInfo *s, Pdata_pos *p, Pdata_pos *moved, 
              float *si, double boxsize, long long ngrid, int s_intp)  {
  //-> move particles by 'si', stored in 'moved'<<- //
  long long i, ip, m, n, l, m1, n1, l1, idx[3], npart, pid;
  float moved_pos, dsi, dp[3], v[2][2][2], xmin, xmax, dx;

  npart=ngrid*ngrid*ngrid;

  // ->> box boundary <<- //
  xmin=0.; xmax=boxsize;
  dx=(xmax-xmin)/(float)ngrid;
       
  printf("->> moving particles ... \n");

  #ifdef _OMP_
  #pragma omp parallel for private(ip,i,m,n,l,m1,n1,l1,idx,moved_pos,dsi,dp,v)
  #endif
  for(ip=0; ip<npart; ip++) {

    // ->> position index <<- //
    pid=p[ip].pid-1;
    pidtogrid(pid, ngrid, idx);

    // ->> do not interpolate <<- //
    if(s_intp==FALSE){

      //->> particle-moving DO NOT interpolation.\n");
      for(i=0; i<3; i++){
        moved_pos=p[ip].pos[i]+ArrayAccess4D_n4(si, 3, ngrid, ngrid, ngrid, i, idx[0], idx[1], idx[2]);

        // ->> periodic boundary condition <<- //
	if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
	else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
	else{
	  moved[ip].pos[i]=moved_pos; }
	}
      }

    // ->> interpolate shift field onto particle position <<- //
    else if(s_intp==TRUE) {
      // ->> ("particle-moving interpolation.\n") <<- //
      // ->> distance to grid <<- //
      for(i=0; i<3; i++) {
        dp[i]=(float)((double)p[ip].pos[i]-s->pmin[i]-idx[i]*s->dpart[i]); 
	if(dp[i]>=(float)(s->pmax[i]-s->pmin[i])) dp[i]=0.; 

	if((dp[i]<0)||(dp[i]>(float)s->dpart[i]))  {
	  printf("dp error, dp=%f (pos=%f, pos_grid=%f, idx=%d, dpart=%f)\n", dp[i], p[ip].pos[i], 
	         idx[i]*(float)s->dpart[i], idx[i], (float)s->dpart[i]);  
	  fflush(stdout);
	  }
	}

      for(i=0; i<3; i++)  {

        // ->> define vertices for 3D interpolation <<- //
        for(l=0; l<2; l++) {
	  l1=l+idx[0]; 
	  if(l1>=s->ngrid)  l1=l1-s->ngrid;

          for(m=0; m<2; m++) {
	    m1=m+idx[1]; 
	    if(m1>=s->ngrid)  m1=m1-s->ngrid;

            for(n=0; n<2; n++) {
	      n1=n+idx[2]; 
	      if(n1>=s->ngrid)  n1=n1-s->ngrid;


              v[l][m][n]=ArrayAccess4D_n4(si, 3, s->ngrid, s->ngrid, s->ngrid, i, l1, m1, n1);
	      }
	    }
	  }

        // ->> tri-linear interpolation <<- //
        dsi=trilinear(p[ip].pos, dp, v);
        moved_pos=p[ip].pos[i]+dsi;

	if(moved_pos<0){
          moved[ip].pos[i]=moved_pos+xmax; }
	else if(moved_pos>=xmax){
          moved[ip].pos[i]=moved_pos-xmax; }
	else{
	  moved[ip].pos[i]=moved_pos; }
        }
      }

    else{abort();}

    }


  //#define _MOVED_PART_OUTPUT_
  #ifdef _MOVED_PART_OUTPUT_
  FILE *fp=fopen(s->test_fname, "wb");

  for(ip=0; ip<s->npart; ip++) 
    fwrite(moved[ip].pos, sizeof(float), 3, fp);
 
  fwrite(si, sizeof(float), s->ngrid*s->ngrid*s->ngrid*3, fp);

  fclose(fp);
  #endif


  printf("->> particles displacement is done.\n");
  return;
  }



void move_grid_general(SimInfo *s, Pdata_pos *moved, float *si) {
  //->> grid moving, no need to generate the grid <<- //
  long long i, j, k, m, ip, grid_min, grid_max;
  float grid[3], xmin, xmax, dx, moved_pos;

  // ->> box boundary <<- //
  xmin=0.;//(s->boxsize-s->bsize_trim)/2.; 
  xmax=xmin+s->bsize_trim;
  dx=(xmax-xmin)/(float)s->ng_trim;

  printf("\n->> shifting uniform grid ...\n");

  #ifdef _OMP_
  #pragma omp parallel for private(i,j,k,m,ip,grid,moved_pos)
  #endif
  for(i=0; i<s->ng_trim; i++)
    for(j=0; j<s->ng_trim; j++)
      for(k=0; k<s->ng_trim; k++) {

        // ->> grid index <<- //
        ip=MemIdx3D(s->ng_trim, i, j, k);

        grid[0]=xmin+i*dx;
        grid[1]=xmin+j*dx;
        grid[2]=xmin+k*dx;

        for(m=0; m<3; m++)  {
          moved_pos=grid[m]+ArrayAccess2D_n2(si, 3, s->npart_trim, m, ip);

          // ->> periodic boundary condition <<- //
          if(moved_pos<0) {
            moved[ip].pos[m]=moved_pos+xmax; }
          else if(moved_pos>=xmax){
            moved[ip].pos[m]=moved_pos-xmax; }
          else{
            moved[ip].pos[m]=moved_pos; }
          }
      }

  printf("->> shifting grid is done.\n");
  return;
  }




