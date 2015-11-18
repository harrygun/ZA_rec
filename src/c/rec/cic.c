#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./nrsrc/nrutil.h"


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
 
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];
} header1;

long long    NumPart, Ngas, NDM;

struct particle_data 
{
  float  Pos[3];
} *P;

float Pos[3],Vel[3];

double  Time, Redshift;

float ***delta;

long long Ngridx,Ngridy,Ngridz;

float xmin,xmax,ymin,ymax,zmin,zmax;

float rhoc,rhom;





/* here the particle data is at your disposal 
 */
int do_what_you_want(void)
{
   long long ip;
   float dx,dy,dz;
   FILE *fp;
   int i,j,k;
   char fname[200];

   dx=(xmax-xmin)/Ngridx;
   dy=(ymax-ymin)/Ngridy;
   dz=(zmax-zmin)/Ngridz;

   printf("dx %f dy %f dz %f\n",dx,dy,dz);

   for(ip=1;ip <= NDM;ip++)
     {
      P[ip].Pos[0]=(P[ip].Pos[0]-xmin)/dx;
      P[ip].Pos[1]=(P[ip].Pos[1]-ymin)/dy;
      P[ip].Pos[2]=(P[ip].Pos[2]-zmin)/dz;
     }

   printf("Ngridx %ld Ngridy %ld Ngridz %ld\n",Ngridx,Ngridy,Ngridz);
   
   printf("density...\n");
   density();
/* write a binary file here */

   sprintf(fname,"%s","delta.binary");
   printf("%s%s\n","writing ...",fname);
   fp=fopen(fname,"wb");
    for(k=0;k<=Ngridx-1;k++)
       for(j=0;j<=Ngridy-1;j++)
          for(i=0;i<=Ngridz-1;i++)
             fwrite(&delta[i][j][k],sizeof(int),1,fp);

   fclose(fp);
   printf("binary done\n"); 
}

void density(void)
{
        long long ip;
	long long i,j,k,i1,j1,k1;
        float xc,yc,zc,dx,dy,dz,tx,ty,tz,x1,y1,z1;
        double aa,masstot;

        for (i=0;i<=Ngridx-1;i++)
            for(j=0;j<=Ngridy-1;j++)
               for(k=0;k<=Ngridz-1;k++)
                  delta[i][j][k]=0.0;

        masstot=0.0;

        for (ip=1;ip<=NDM;ip++) 
            {
            i=(int)P[ip].Pos[0]; xc=(float)i;
            j=(int)P[ip].Pos[1]; yc=(float)j;
            k=(int)P[ip].Pos[2]; zc=(float)k;

            if(i<0) i=i+Ngridx;
            if(j<0) j=j+Ngridy;
            if(k<0) k=k+Ngridz;
      
            if(i>=Ngridx) i=i-Ngridx;
            if(j>=Ngridy) j=j-Ngridy;
            if(k>=Ngridz) k=k-Ngridz;

            dx=fabs(P[ip].Pos[0]-xc); tx=fabs(1.0-dx);
            dy=fabs(P[ip].Pos[1]-yc); ty=fabs(1.0-dy);
            dz=fabs(P[ip].Pos[2]-zc); tz=fabs(1.0-dz);

            i1=i+1;  
            j1=j+1; 
	    k1=k+1;

            if(i1<0) i1=i1+Ngridx;
            if(j1<0) j1=j1+Ngridy;
            if(k1<0) k1=k1+Ngridz;

            if(i1>=Ngridx) i1=i1-Ngridx;
            if(j1>=Ngridy) j1=j1-Ngridy;
            if(k1>=Ngridz) k1=k1-Ngridz;

            delta[i][j][k]+=header1.mass[1]*tx*ty*tz;
            delta[i1][j][k]+=header1.mass[1]*dx*ty*tz;
            delta[i][j1][k]+=header1.mass[1]*tx*dy*tz;
            delta[i1][j1][k]+=header1.mass[1]*dx*dy*tz;
            delta[i][j][k1]+=header1.mass[1]*tx*ty*dz;
            delta[i1][j][k1]+=header1.mass[1]*dx*ty*dz;
            delta[i][j1][k1]+=header1.mass[1]*tx*dy*dz;
            delta[i1][j1][k1]+=header1.mass[1]*dx*dy*dz;
            masstot+=header1.mass[1]*tx*ty*tz;
	    masstot+=header1.mass[1]*dx*ty*tz;
	    masstot+=header1.mass[1]*tx*dy*tz;
	    masstot+=header1.mass[1]*dx*dy*tz;
	    masstot+=header1.mass[1]*tx*ty*dz;
	    masstot+=header1.mass[1]*dx*ty*dz;
	    masstot+=header1.mass[1]*tx*dy*dz;
	    masstot+=header1.mass[1]*dx*dy*dz;
	}	

      printf("masstot summed: %lg masstot all: %lg\n",masstot*1e10/header1.HubbleParam,header1.mass[1]*1e10/header1.HubbleParam*NDM);

       aa=0;
       x1=(xmax-xmin)/Ngridx/header1.HubbleParam;
       y1=(ymax-ymin)/Ngridy/header1.HubbleParam;
       z1=(zmax-zmin)/Ngridz/header1.HubbleParam;
       for(i=0;i<=Ngridx-1;i++)
          for(j=0;j<=Ngridy-1;j++)
             for(k=0;k<=Ngridz-1;k++){
                delta[i][j][k]=delta[i][j][k]/(x1*y1*z1)/rhom*1e10/header1.HubbleParam-1;
                aa+=delta[i][j][k];
                }
       printf("mean density %lg\n",aa/Ngridx/Ngridy/Ngridz);
}







