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

  #include "glbvarb.h"

  #include "cic.h"


#ifdef _MPI_
  #include "mpi.h"
#endif

#ifdef _OPEN_MP_
  #include "omp.h"
#endif


  int main( int argc, char *argv[])  {

      int debug= 20, i, j, do_MPI=99;
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


      char *smooth_type, *prop_file_name, *output_prefix, *plin_name;
      int  do_propagator;

      smooth_type = iniparser_getstring(dict, "LPT_LOG:smooth_type", "gaussian");
      cp.R = iniparser_getdouble(dict, "LPT_LOG:smooth_scale", 3.) ;
      cp.z = iniparser_getdouble(dict, "LPT_LOG:redshift", 0) ;

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

      do_propagator=iniparser_getboolean
                      (dict, "EPT_LOG_RESUM:do_propagator", INIFALSE);
      prop_file_name=iniparser_getstring
                    (dict,"EPT_LOG_RESUM:propagator_file_name", "fkchart.dat");
      output_prefix=iniparser_getstring
                    (dict,"EPT_LOG_RESUM:output_prefix", NULL);

      plin_name =iniparser_getstring
                    (dict,"EPT_LOG_RESUM:linear_power_name", NULL);

      /*-----------------------------------------------------------------------*/

      Interpar *power = (Interpar *)malloc(sizeof(Interpar));
      FILE *fp = fopen(plin_name, "r"); 

      init_powerInterp(power, fp);
      cp.power=power;

      cp.zinit = iniparser_getdouble(dict, "EPT_LOG_RESUM:initial_redshift", 35);
      cp.a = Dp(&cp, cp.z) / Dp(&cp, cp.zinit);
      printf("zinit=%lg\n", cp.zinit);

  /*-----------------------------------------------
              End of initialization.
  -----------------------------------------------*/


  /*--------------------------
    Initialize the data
  ----------------------------*/







/*-----------------------------------------------------
                   free all
-----------------------------------------------------*/
      iniparser_freedict(dict);


  #ifdef _MPI_
//      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
  #endif

      
    }


