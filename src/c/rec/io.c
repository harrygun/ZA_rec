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

int snap_number;

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */

/* modified from read_snapshot.c in the Gadget package*/

int main(int argc, char **argv)
{
  char path[200], input_fname[200],basename[200];
  int  snapshot_number, files;
  int i,j,k;

  Ngridx=100;
  Ngridy=100;
  Ngridz=100;

  delta=f3tensor(0L,Ngridx-1,0L,Ngridy-1,0L,Ngridz-1);
  for(i=0;i<Ngridx;i++)
     for(j=0;j<Ngridy;j++)
        for(k=0;k<Ngridz;k++)
           delta[i][j][k]=0.0;

  files=1;                             

  //sprintf(input_fname,"%s","../ics_L50_N100_0_0_00.0.data");

  sprintf(input_fname,"%s","/home1/yuebin/Data/FirstStars/coarse512/init_gadget_c.dat");
  printf("Ngridx %ld Ngridy %ld Ngridz %ld\n",Ngridx,Ngridy,Ngridz);

  load_snapshot(input_fname, files);

  do_what_you_want();
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  long long    t,n,off,pc,pc_new,pc_sph,ip;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
 
 for(i=0, pc=1; i<files; i++, pc=pc_new)
     {
	if(files>1)
	    sprintf(buf,"%s.%d",fname,i);
	else
	    sprintf(buf,"%s",fname);		
	
        if(!(fd=fopen(buf,"r")))
	   {
	    printf("can't open file `%s`\n",buf);
	    exit(0);
	   }
		
        printf("reading `%s' ...\n",buf); fflush(stdout);
		
	fread(&dummy, sizeof(dummy), 1, fd);
	fread(&header1, sizeof(header1), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);
	
	for(k=0;k<6;k++)
	   {
printf("Number of Type %d Particles: %d, mass %lf [/h M_sun]\n",k,header1.npart[k],header1.mass[k]*1e10);
	   }    
	printf("BoxSize:%f\n",header1.BoxSize);
	printf("Redshift:%f\n",header1.redshift);
	 

	for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	   {
	    NumPart+=(unsigned int) header1.npartTotal[k];
	    NumPart+=(((long long) header1.npartTotalHighWord[k]) << 32);
	   }

       NDM=NumPart;
       if(i==0)
         allocate_memory(); 
		
        printf("NumPart %ld NDM %ld\n",NumPart,NDM);
		
		
	SKIP;
	for(k=0,pc_new=pc;k<6;k++)
	   {
	   for(n=0;n<header1.npart[k];n++)
	      {
	      fread(&Pos[0], sizeof(float), 3, fd);
              if(k==1)
	        {			
		P[pc_new].Pos[0]=Pos[0];
	        P[pc_new].Pos[1]=Pos[1];
		P[pc_new].Pos[2]=Pos[2];
		pc_new++;
	        }
	      }
	}

        SKIP;
	fclose(fd);
		
    }
	
  printf("read done\n");

  Time= header1.time;
  Redshift= header1.time;
  rhoc=(2.7752e11)*header1.HubbleParam*header1.HubbleParam;
  rhom=rhoc*header1.Omega0;

  xmin=header1.BoxSize/2.0;
  xmax=header1.BoxSize/2.0;

  ymin=header1.BoxSize/2.0;
  ymax=header1.BoxSize/2.0;

  zmin=header1.BoxSize/2.0;
  zmax=header1.BoxSize/2.0;

  for(ip=1;ip <= NDM;ip++)
     {
       if(P[ip].Pos[0] < xmin) xmin=P[ip].Pos[0];
       if(P[ip].Pos[0] > xmax) xmax=P[ip].Pos[0];

       if(P[ip].Pos[1] < ymin) ymin=P[ip].Pos[1];
       if(P[ip].Pos[1] > ymax) ymax=P[ip].Pos[1];

       if(P[ip].Pos[2] < zmin) zmin=P[ip].Pos[2];
       if(P[ip].Pos[2] > zmax) zmax=P[ip].Pos[2];
     }
 printf("xmin %f ymin %f zmin %f\n",xmin,ymin,zmin);
 printf("xmax %f ymax %f zmax %f\n",xmax,ymax,zmax);
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");
  printf("NDM %ld\n",NDM);

  if(!(P=malloc(NDM*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */
  printf("allocating memory...done\n");
}




  











