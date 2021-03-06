  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>
  #include <omp.h>

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
  #include "poisson.h"
  #include "reconstruction_partmoving.h"



#ifdef _MPI_
  #include "mpi.h"
#endif






void testing_fftw(SimInfo *s, float *d, float *phi, float *phi_i, float *phi_ij, int fft_return_type);



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

      // simulation info & reconstruction constroller //
      SimInfo s;
      RectCtrl rc;

      //double boxsize;
      char *smooth_type, *particle_fname, *droot, *plin_name;
      char *oden_fname, *main_dtype, *rec_name;
      int do_density, save_odensity;

      cp.R = iniparser_getdouble(dict, "Rect:smooth_R", 10.);
      cp.z = iniparser_getdouble(dict, "Rect:redshift", 0) ;

      // ->> 
      s.smooth_R=cp.R;
      s.smooth_type = iniparser_getstring(dict, "Rect:smooth_type", "Gaussian");
      s.boxsize = iniparser_getdouble(dict, "Rect:boxsize", 0) ;
      s.ngrid=iniparser_getint(dict, "Rect:ngrid", 0);
      s.npart=pow(s.ngrid,3);
      printf("boxsize=%lg, ngrid=%d, npart=%d\n",s.boxsize,s.ngrid,s.npart);

      if(strcmp( s.smooth_type, "Tophat")==0 ) {
        //cp.flg[0]= _TOPHAT_SMOOTH_ ;
        s.smooth_type_flag= _TOPHAT_SMOOTH_ ;
	printf("Tophat smoothing.\n");
	}
      else if(strcmp( s.smooth_type, "Gaussian")==0 ) {
        //cp.flg[0]= _GAUSSIAN_SMOOTH_ ;
        s.smooth_type_flag= _GAUSSIAN_SMOOTH_ ;
	printf("Gussian smoothing.\n");
	}
      else
        myerr("smooth type error", FALSE);

      printf("smooth scale=%lg\n", s.smooth_R);
      printf("z=%lg\n", cp.z);

      droot=iniparser_getstring(dict,"Rect:data_root", "~/");
      particle_fname=iniparser_getstring(dict,"Rect:particle_file_name", "x.dat");
      //printf("particle fname: %s\n", particle_fname);

      /*-----------------------------------------------------------------------*/
      // ->> read controller <<- // 
      rc.do_rect=iniparser_getboolean(dict, "Rect:do_reconstruction", INIFALSE);
      rc.displacement_intp=iniparser_getboolean(dict, "Rect:displacement_interpolation", INIFALSE);
      rc.do_disp_perturb=iniparser_getboolean(dict, "Rect:perturbe_displacement", INIFALSE);
      //->> reconstruction type, displacement type etc. <<- //
      rc.rec_type=iniparser_getstring(dict, "Rect:reconstruction_type", "displaced");
      rc.displacement_type=iniparser_getstring(dict,"Rect:displacement_type", "backward");
      rc.displacement_order=iniparser_getstring(dict,"Rect:displacement_order", "1LPT");

      // ->> other controller <<- //
      do_density=iniparser_getboolean(dict, "Rect:do_density", INIFALSE);
      save_odensity=iniparser_getboolean(dict, "Rect:save_original_density", INIFALSE);
      oden_fname=iniparser_getstring(dict,"Rect:original_density_fname", "y.dat");
      rec_name=iniparser_getstring(dict,"Rect:reconstructed_fname", "y.dat");
      //->>
      snprintf(rc.rec_fname, 200, "%s_%s_R%g.dat", rec_name, s.smooth_type, s.smooth_R);

      /*-----------------------------------------------------------------------*/
      // ->> initialize power spectrum <<- //
      plin_name=iniparser_getstring(dict,"Rect:linear_power_name", NULL);
      Interpar *power = (Interpar *)malloc(sizeof(Interpar));
      FILE *fp = fopen(plin_name, "r"); 
      init_powerInterp(power, fp);
      cp.power=power;
      fclose(fp);

      /*-----------------------------------------------------------------------*/
      int do_fftw_testing;
      char *fftw_test_fname, *fftw_return_type;
      do_fftw_testing=iniparser_getboolean(dict, "Rect:do_fftw_testing", INIFALSE);
      fftw_test_fname=iniparser_getstring(dict,"Rect:fftw_test_fname", "y.dat");
      fftw_return_type=iniparser_getstring(dict,"Rect:fftw_test_return_type", "gradient");

      s.test_fname=iniparser_getstring(dict,"Rect:other_test_fname", "y.dat");
    /*-----     End of initialization.    ------*/


    /*-----------------------------------------------------
            // ->>     Importing data     <<- //
    -----------------------------------------------------*/
    //  ->> loading particle data <<- //
    Pdata_pos *p=(Pdata_pos *)malloc(s.npart*sizeof(Pdata));
    load_cita_simulation_position(particle_fname, p, s.npart);

    float *d;
    double rhom_, dmean; //, particle_mass;

    for(i=0; i<3; i++){ s.ngrid_xyz[i]=s.ngrid; }
    d=(float *)malloc(sizeof(float)*s.ngrid_xyz[0]*s.ngrid_xyz[1]*s.ngrid_xyz[2]);

    // ->> particle mass <<- //
    rhom_=get_rhom(&cp, cp.z);
    s.particle_mass=part_mass(&cp, cp.z, s.boxsize, s.ngrid);
    printf("\nmean density=%lg (M_star)*(h/Mpc)^3\nparticle mass resolution=%lg\n\n", rhom_, s.particle_mass);

    // ->> if do CIC density estimation <<- //
    if (do_density==TRUE) {
      // ->> CIC density estimation <<- //
      dmean=cic_density(p, d, s.boxsize, s.particle_mass, s.npart, s.ngrid_xyz, NULL); 

      // ->> save density field in file <<- //
      if(save_odensity==TRUE) { //
        fp=fopen(oden_fname, "wb");
        for(i=0; i<s.ngrid; i++)
          for(j=0; j<s.ngrid; j++)
            for(k=0; k<s.ngrid; k++)
              fwrite(&ArrayAccess3D(d, s.ngrid, i, j, k), sizeof(float), 1, fp);
        fclose(fp);
        }

      }
    else { 
      // ->> otherwise, import density field <<- //
      main_dtype="float";
      printf("Importing the original density map directly, dtype=%s\n", main_dtype);
      load_scalar_map(oden_fname, d, s.ngrid, main_dtype);
      }

    if(do_fftw_testing==TRUE) {
      printf("Do FFTW testing.\n"); fflush(stdout);
      fftw_tester(&s, d, fftw_return_type, fftw_test_fname);
      abort();
      }


    /*-----------------------------------------------------
         // ->>   performing reconstruction   <<- //
    -----------------------------------------------------*/
    if(rc.do_rect!=TRUE) {
      printf("Do NOT perform reconstruction.\n");
      fflush(stdout); abort(); }


    // ->> Perform some test <<- //
    int fft_return_type;
    float *phi, *phi_i, *phi_ij;
   
    //phi=(float *)fftwf_malloc(sizeof(float)*pow(s.ngrid,3));
    phi_i=(float *)fftwf_malloc(sizeof(float)*pow(s.ngrid,3)*3);
    //phi_ij=(float *)fftwf_malloc(sizeof(float)*pow(s.ngrid,3)*9);
          
    
    //fft_return_type=_RETURN_PHI_ONLY_;
    fft_return_type=_RETURN_GRADIENT_;
    testing_fftw(&s, d, phi, phi_i, phi_ij, fft_return_type);


    // ->> write files <<- //
    fp=fopen(s.test_fname, "wb");
    fwrite(phi_i, sizeof(float), s.ngrid*s.ngrid*s.ngrid*3, fp);

    fclose(fp);




    /*-----------------------------------------------------
                       free all
    -----------------------------------------------------*/
    iniparser_freedict(dict);
    free(power);
    free(p); free(d);

    //fftwf_free(phi);

    fftwf_free(phi_i);
    //fftwf_free(phi_ij);



    #ifdef _MPI_
      MPI_Finalize();
    #endif

    }





void testing_fftw(SimInfo *s, float *d, float *phi, float *phi_i, float *phi_ij, int fft_return_type){
  
  poisson_solver_float(d, phi, phi_i, phi_ij, s->boxsize, s->ngrid, s->smooth_type_flag, s->smooth_R, fft_return_type);

  return;
  }
