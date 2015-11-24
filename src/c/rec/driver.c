  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <iniparser.h>

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
  #include "cic.h"


#ifdef _MPI_
  #include "mpi.h"
#endif

#ifdef _OPEN_MP_
  #include "omp.h"
#endif


  int main( int argc, char *argv[])  {

      int debug= 20, i, j, k, do_MPI=99;
      char *ini_name;

      if(argc!=2)  //
        myerr("Input parameters file is needed.", FALSE);

        ini_name=argv[1];


/*------------------------------------------------
            MPI initialization.
------------------------------------------------*/
    #ifdef _MPI_
      do_MPI = TRUE;

      int mpi_ntask, mpi_rank, mpi_rc;
      mpi_rc = MPI_Init(&argc, &argv);

      if (mpi_rc != MPI_SUCCESS) {
           printf ("Error starting MPI program. Terminating.\n");
           MPI_Abort(MPI_COMM_WORLD, mpi_rc);
           }

      MPI_Comm_size(MPI_COMM_WORLD,&mpi_ntask);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

      printf ("Number of tasks= %d My rank= %d\n", mpi_ntask, mpi_rank);

  /* MPI initialization ends. */
      if(mpi_rank==0) 
        printf("%d Sending parameter filename %s to other processes.\n", mpi_rank, ini_name);

      MPI_Bcast( ini_name, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

      if(mpi_rank!=0) 
        printf("%d Received parameter filename %s.\n", mpi_rank, ini_name);

    #endif

    #ifndef _MPI_
      do_MPI = FALSE;
    #endif

  /*--------------------------------------
         End of MPI initialization.
  --------------------------------------*/



  /*--------------------------------------------------*/
      printf("Opening File '%s'.\n", ini_name);

      dictionary * dict;
      dict = iniparser_load(ini_name);
    /*-----------------------------------------------
          initialize the cosmological parameters
    -----------------------------------------------*/

      Cospar  cp;
      init_cospar(&cp, dict);

      if(debug>=20) {
        printf("omem=%lg\n", cp.omem);
        printf("omek=%lg\n", cp.omek);
        printf("omeb=%lg\n", cp.omeb);
        printf("h0=%lg\n", cp.h0);
        printf("w0=%lg\n", cp.w0);
        printf("w1=%lg\n", cp.w1);
        }

      double boxsize;
      char *smooth_type, *particle_fname, *droot, *plin_name, *oden_fname;
      int do_density, do_potential, do_rect, save_odensity;
      int npart, ngrid;

      smooth_type = iniparser_getstring(dict, "Rect:smooth_type", "gaussian");
      cp.R = iniparser_getdouble(dict, "Rect:smooth_scale", 10.);
      cp.z = iniparser_getdouble(dict, "Rect:redshift", 0) ;


      boxsize = iniparser_getdouble(dict, "Rect:boxsize", 0) ;
      ngrid=iniparser_getint(dict, "Rect:grid_size", 0);
      npart=pow(ngrid,3);
      printf("boxsize=%lg,  ngrid=%d,  npart=%d\n", boxsize, ngrid, npart);

      if(strcmp( smooth_type, "tophat")==0 ) {
        cp.flg[0]= _TOPHAT_SMOOTH_ ;
	printf("Tophat smoothing.\n");
	}
      else if(strcmp( smooth_type, "gaussian")==0 ) {
        cp.flg[0]= _GAUSSIAN_SMOOTH_ ;
	printf("Gussian smoothing.\n");
	}
      else
        myerr("smooth type error", FALSE);

      printf("smooth scale=%lg\n", cp.R);
      printf("z=%lg\n", cp.z);

      droot=iniparser_getstring(dict,"Rect:data_root", "~/");
      particle_fname=iniparser_getstring(dict,"Rect:particle_file_name", "x.dat");
      //printf("particle fname: %s\n", particle_fname);

      /*-----------------------------------------------------------------------*/
      // ->> read controller <<- // 
      do_rect=iniparser_getboolean(dict, "Rect:do_reconstruction", INIFALSE);
      do_density=iniparser_getboolean(dict, "Rect:do_density", INIFALSE);
      do_potential=iniparser_getboolean(dict, "Rect:do_potential", INIFALSE);

      save_odensity=iniparser_getboolean(dict, "Rect:save_original_density", INIFALSE);
      oden_fname=iniparser_getstring(dict,"Rect:original_density_fname", "y.dat");

      /*-----------------------------------------------------------------------*/
      // ->> initialize power spectrum <<- //
      plin_name=iniparser_getstring(dict,"Rect:linear_power_name", NULL);
      Interpar *power = (Interpar *)malloc(sizeof(Interpar));
      FILE *fp = fopen(plin_name, "r"); 
      init_powerInterp(power, fp);
      cp.power=power;
      fclose(fp);

    /*-----     End of initialization.    ------*/



    /*-----------------------------------------------------
            // ->>     Importing data     <<- //
    -----------------------------------------------------*/
    //  ->> loading particle data <<- //
    Pdata_pos *p=(Pdata_pos *)malloc(npart*sizeof(Pdata));
    load_cita_simulation_position(particle_fname, p, npart);


   // ->> if do CIC density estimation <<- //
   float ***d;
   double particle_mass, rhom_, dmean;
   int ngrid_xyz[3];
    if (do_density=TRUE) {
      // ->> 
      for(i=0; i<3; i++)
        ngrid_xyz[i]=ngrid;

      d=(float ***)anymat3(ngrid_xyz[0], ngrid_xyz[1], ngrid_xyz[2], 
                           sizeof(float),sizeof(float *),sizeof(float **));

      rhom_=get_rhom(&cp, cp.z);
      particle_mass=part_mass(&cp, cp.z, boxsize, ngrid);
      printf("\nmean density=%lg (M_star)*(h/Mpc)^3, particle mass resolution=%lg\n\n", rhom_, particle_mass);

      // ->> CIC density estimation <<- //
      dmean=cic_density(p, d, boxsize, particle_mass, npart, ngrid_xyz); 

      // ->> save density field in file <<- //
      if(save_odensity==TRUE) { //
        fp=fopen(oden_fname, "wb");
        for(i=0; i<ngrid; i++)
          for(j=0; j<ngrid; j++)
            for(k=0; k<ngrid; k++)
              fwrite(&d[i][j][k], sizeof(float), 1, fp);
        fclose(fp);
        }

      }
    // ->> otherwise, import density field <<- //
    else { }

    /*-----------------------------------------------------
         // ->>   performing reconstruction   <<- //
    -----------------------------------------------------*/

    // ->> Obtain displacement field <<- //




    /*-----------------------------------------------------
                       free all
    -----------------------------------------------------*/
    iniparser_freedict(dict);
    free(p);
    freemat3((void ***)d, ngrid_xyz[0], ngrid_xyz[1], ngrid_xyz[2]);


    #ifdef _MPI_
      MPI_Finalize();
      //MPI_Barrier(MPI_COMM_WORLD);
    #endif

    }


