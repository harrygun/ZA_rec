  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>


  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "parvar.h"
  #include "io.h"

  #ifdef _OMP_
  #include <omp.h>
  #endif



double get_rhom(Cospar *cp, double z){
  //->>  M_star*(h/Mpc)^3 <<-//
  return (2.7752e11)*cp->omem/cp->h0*pow(1.+z,3.);
  }

double part_mass(Cospar *cp, double z, double boxsize, int ngrid){
  double rhom=get_rhom(cp, z);
  return rhom*pow(boxsize, 3.)/pow((double)(ngrid), 3.);
  }



double density(Pdata_pos *p, float *d, double mass, double pmin[3], double pmax[3], double dg[3],
             long long npart, long long ngridx, long long ngridy, long long ngridz) {
  long long ip, i,j,k,i1,j1,k1, n;
  double xc,yc,zc,dx,dy,dz,tx,ty,tz,x1,y1,z1;
  double aa, masstot, dmean, pos_norm[3];
  
  // ->> 
  for(i=0; i<ngridx; i++)
    for(j=0; j<ngridy; j++)
      for(k=0; k<ngridz; k++)
        //d[i][j][k]=0.0;
        ArrayAccess3D(d, ngridx, i, j, k)=0.0;
  masstot=0.0;

  #ifdef _OMP_
  #pragma omp parallel for private(ip,i,j,k,i1,j1,k1,n,xc,yc,zc,pos_norm,tx,ty,tz,dx,dy,dz)
  #endif
  for(ip=0; ip<npart; ip++) {

    // ->> renormalize 
    for(n=0; n<3; n++)
      pos_norm[n]=(p[ip].pos[n]-pmin[n])/dg[n];

    //i=(int)p[ip].pos[0]; xc=(double)i;
    //j=(int)p[ip].pos[1]; yc=(double)j;
    //k=(int)p[ip].pos[2]; zc=(double)k;

    i=(long long)pos_norm[0]; xc=(double)i;
    j=(long long)pos_norm[1]; yc=(double)j;
    k=(long long)pos_norm[2]; zc=(double)k;
    
    if(i<0) i=i+ngridx;
    if(j<0) j=j+ngridy;
    if(k<0) k=k+ngridz;
    
    if(i>=ngridx) i=i-ngridx;
    if(j>=ngridy) j=j-ngridy;
    if(k>=ngridz) k=k-ngridz;
    
    //dx=fabs(p[ip].pos[0]-xc); tx=fabs(1.0-dx);
    //dy=fabs(p[ip].pos[1]-yc); ty=fabs(1.0-dy);
    //dz=fabs(p[ip].pos[2]-zc); tz=fabs(1.0-dz);


    dx=fabs(pos_norm[0]-xc); tx=fabs(1.0-dx);
    dy=fabs(pos_norm[1]-yc); ty=fabs(1.0-dy);
    dz=fabs(pos_norm[2]-zc); tz=fabs(1.0-dz);
    
    i1=i+1;  
    j1=j+1; 
    k1=k+1;
    
    if(i1<0) i1=i1+ngridx;
    if(j1<0) j1=j1+ngridy;
    if(k1<0) k1=k1+ngridz;
    
    if(i1>=ngridx) i1=i1-ngridx;
    if(j1>=ngridy) j1=j1-ngridy;
    if(k1>=ngridz) k1=k1-ngridz;

    /*
    d[i][j][k]+=mass*tx*ty*tz;
    d[i1][j][k]+=mass*dx*ty*tz;
    d[i][j1][k]+=mass*tx*dy*tz;
    d[i1][j1][k]+=mass*dx*dy*tz;
    d[i][j][k1]+=mass*tx*ty*dz;
    d[i1][j][k1]+=mass*dx*ty*dz;
    d[i][j1][k1]+=mass*tx*dy*dz;
    d[i1][j1][k1]+=mass*dx*dy*dz;
    */

    ArrayAccess3D(d, ngridx, i, j, k) +=mass*tx*ty*tz;
    ArrayAccess3D(d, ngridx, i1, j, k)+=mass*dx*ty*tz;
    ArrayAccess3D(d, ngridx, i, j1, k) +=mass*tx*dy*tz;
    ArrayAccess3D(d, ngridx, i1, j1, k) +=mass*dx*dy*tz;
    ArrayAccess3D(d, ngridx, i, j, k1) +=mass*tx*ty*dz;
    ArrayAccess3D(d, ngridx, i1, j, k1) +=mass*dx*ty*dz;
    ArrayAccess3D(d, ngridx, i, j1, k1) +=mass*tx*dy*dz;
    ArrayAccess3D(d, ngridx, i1, j1, k1)+=mass*dx*dy*dz;


    masstot+=mass*tx*ty*tz;
    masstot+=mass*dx*ty*tz;
    masstot+=mass*tx*dy*tz;
    masstot+=mass*dx*dy*tz;
    masstot+=mass*tx*ty*dz;
    masstot+=mass*dx*ty*dz;
    masstot+=mass*tx*dy*dz;
    masstot+=mass*dx*dy*dz;
  }

  aa=0;
  x1=(pmax[0]-pmin[0])/(double)ngridx;
  y1=(pmax[1]-pmin[1])/(double)ngridy;
  z1=(pmax[2]-pmin[2])/(double)ngridz;

  for(i=0; i<ngridx; i++)
    for(j=0; j<ngridy; j++)
      for(k=0; k<ngridz; k++) {
        //d[i][j][k]=d[i][j][k]/(x1*y1*z1);
        //aa+=d[i][j][k];

        ArrayAccess3D(d, ngridx, i, j, k)=ArrayAccess3D(d, ngridx, i, j, k)/(x1*y1*z1);
        aa+=ArrayAccess3D(d, ngridx, i, j, k);
        }
	
  dmean=aa/(double)ngridx/(double)ngridy/(double)ngridz;
  printf("mean density %lg\n", dmean);

  //->>
  for(i=0; i<ngridx; i++)
    for(j=0; j<ngridy; j++)
      for(k=0; k<ngridz; k++) {
        //d[i][j][k]=d[i][j][k]/dmean-1.;
        ArrayAccess3D(d, ngridx, i, j, k)=ArrayAccess3D(d, ngridx, i, j, k)/dmean-1.;
        }

  return dmean;
  }





double cic_density(Pdata_pos *p, float *d, double boxsize, 
                      double mass, long long npart, long long ngrid[3]) {
  long long ip;
  double xmin, ymin, zmin, xmax, ymax, zmax, pmin[3], pmax[3], dx[3], dmean;

  // ->> obtain boundary <<- //
  xmin=boxsize/2.0; xmax=boxsize/2.0;
  ymin=boxsize/2.0; ymax=boxsize/2.0;
  zmin=boxsize/2.0; zmax=boxsize/2.0;

  for(ip=0; ip<npart; ip++) {
    if(p[ip].pos[0]<xmin) xmin=p[ip].pos[0];
    if(p[ip].pos[0]>xmax) xmax=p[ip].pos[0];
    
    if(p[ip].pos[1]<ymin) ymin=p[ip].pos[1];
    if(p[ip].pos[1]>ymax) ymax=p[ip].pos[1];
    
    if(p[ip].pos[2]<zmin) zmin=p[ip].pos[2];
    if(p[ip].pos[2]>zmax) zmax=p[ip].pos[2];
    }

  //printf("xmin %f ymin %f zmin %f\n",xmin, ymin, zmin);
  //printf("xmax %f ymax %f zmax %f\n",xmax, ymax, zmax);


  pmin[0]=xmin; pmin[1]=ymin; pmin[2]=zmin; 
  pmax[0]=xmax; pmax[1]=ymax; pmax[2]=zmax; 


  // ->> dx, dy, dz <<- //
  dx[0]=(xmax-xmin)/(double)ngrid[0];
  dx[1]=(ymax-ymin)/(double)ngrid[1];
  dx[2]=(zmax-zmin)/(double)ngrid[2];
  
  //printf("cic_density: dx=%f, dy=%f, dz=%f\n",dx[0], dx[1], dx[2]);
  
  // ->> renormalize particle positions <<- //
  /*
  for(ip=0; ip<npart; ip++) {
    p[ip].pos[0]=(p[ip].pos[0]-xmin)/dx;
    p[ip].pos[1]=(p[ip].pos[1]-ymin)/dy;
    p[ip].pos[2]=(p[ip].pos[2]-zmin)/dz;
    }
  */
  //printf("ngridx %ld ngridy %ld ngridz %ld\n",ngrid[0],ngrid[1],ngrid[2]);

  // ->> CIC density <<- //
  dmean=density(p, d, mass, pmin, pmax, dx, npart, ngrid[0], ngrid[1], ngrid[2]); 

  printf("\n->> CIC density is done.\n");
  return;
  }
