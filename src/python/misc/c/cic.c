#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//include "./nrsrc/nrutil.h"



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */

/* modified from read_snapshot.c in the Gadget package*/





double density(float *Pos, float *delta, long long NDM, long long Ngrid[3], double mass) 
{

        long long ip, Ngridx, Ngridy, Ngridz;
	long long i,j,k,i1,j1,k1,p_ncol, d_ngrid;
        float xc,yc,zc,dx,dy,dz,tx,ty,tz,x1,y1,z1;
        double aa,masstot;

        Ngridx=Ngrid[0];
        Ngridy=Ngrid[1];
        Ngridz=Ngrid[2];


        masstot=0.0;
	p_ncol=3;

        for (ip=1;ip<=NDM;ip++) 
            {
            i=(int)Pos[ip*ncol+0]; xc=(float)i;
            j=(int)Pos[ip*ncol+1]; yc=(float)j;
            k=(int)Pos[ip*ncol+2]; zc=(float)k;

            if(i<0) i=i+Ngridx;
            if(j<0) j=j+Ngridy;
            if(k<0) k=k+Ngridz;
      
            if(i>=Ngridx) i=i-Ngridx;
            if(j>=Ngridy) j=j-Ngridy;
            if(k>=Ngridz) k=k-Ngridz;

            dx=fabs(Pos[ip*ncol+0]-xc); tx=fabs(1.0-dx);
            dy=fabs(Pos[ip*ncol+1]-yc); ty=fabs(1.0-dy);
            dz=fabs(Pos[ip*ncol+2]-zc); tz=fabs(1.0-dz);

            i1=i+1;  
            j1=j+1; 
	    k1=k+1;

            if(i1<0) i1=i1+Ngridx;
            if(j1<0) j1=j1+Ngridy;
            if(k1<0) k1=k1+Ngridz;

            if(i1>=Ngridx) i1=i1-Ngridx;
            if(j1>=Ngridy) j1=j1-Ngridy;
            if(k1>=Ngridz) k1=k1-Ngridz;

            delta[(i *Ngridy + j )*Ngridz + k ]+=mass*tx*ty*tz;
            delta[(i1*Ngridy + j )*Ngridz + k ]+=mass*dx*ty*tz;
            delta[(i *Ngridy + j1)*Ngridz + k ]+=mass*tx*dy*tz;
            delta[(i1*Ngridy + j1)*Ngridz + k ]+=mass*dx*dy*tz;
            delta[(i *Ngridy + j )*Ngridz + k1]+=mass*tx*ty*dz;
            delta[(i1*Ngridy + j )*Ngridz + k1]+=mass*dx*ty*dz;
            delta[(i *Ngridy + j1)*Ngridz + k1]+=mass*tx*dy*dz;
            delta[(i1*Ngridy + j1)*Ngridz + k1]+=mass*dx*dy*dz;
            masstot+=mass*tx*ty*tz;
	    masstot+=mass*dx*ty*tz;
	    masstot+=mass*tx*dy*tz;
	    masstot+=mass*dx*dy*tz;
	    masstot+=mass*tx*ty*dz;
	    masstot+=mass*dx*ty*dz;
	    masstot+=mass*tx*dy*dz;
	    masstot+=mass*dx*dy*dz;
	}

      //printf("masstot summed: %lg masstot all: %lg\n",masstot*1e10/header1.HubbleParam,header1.mass[1]*1e10/header1.HubbleParam*NDM);


    return masstot;
}







  





