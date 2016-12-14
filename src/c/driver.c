  #include <stdio.h>
  #include <time.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <gsl/gsl_histogram2d.h>
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
  #include "misc.h"
  #include "poisson.h"
  #include "reconstruction_partmoving.h"
  #include "test_likelihood.h"

  #include "stat_model.h"
  #include "maxlike_phi_rec.h"


#ifdef _MPI_
  #include "mpi.h"
#endif

#ifdef _OMP_
  #include <omp.h>
#endif


  int main( int argc, char *argv[])  {

      int debug= 20, i, j, k, do_MPI=99;
      clock_t t0, t1;
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

      if(mpi_rank==0) 
        t0=clock();

    #else
      do_MPI = FALSE;
      t0=clock();
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
      char *smooth_type, *particle_fname, *pid_fname, *droot, *plin_name;
      char *oden_fname, *main_dtype, *rec_name, *fname_pinit, *fname_pid_init, *fname_offset;
      char *raw_disp_fname, *stat_disp_fname, *phi_mlik_fname, *mlik_rec_out_fname;
      int do_density, import_density, save_odensity;

      cp.R = iniparser_getdouble(dict, "Rect:smooth_R", 10.);
      cp.z = iniparser_getdouble(dict, "Rect:redshift", 0) ;
      cp.zinit = iniparser_getdouble(dict, "Rect:initial_redshift", 0) ;

      // ->> 
      s.smooth_R=cp.R;
      s.smooth_type = iniparser_getstring(dict, "Rect:smooth_type", "Gaussian");
      s.boxsize = iniparser_getdouble(dict, "Rect:boxsize", 0) ;
      s.ngrid=iniparser_getint(dict, "Rect:ngrid", 0);
      s.npart=pow(s.ngrid,3);

      // ->> simulation box drift <<- //
      fname_offset=iniparser_getstring(dict, "Rect:simbox_drift_file_name", "None");
      if(strcmp(fname_offset, "None")!=0) {
        load_simulation_offset(fname_offset, s.drift, s.drift_init);
	}
      else{
        s.drift[0]=s.drift[1]=s.drift[2]=0.;
        s.drift_init[0]=s.drift_init[1]=s.drift_init[2]=0.;
        }

      //->> 
      printf("boxsize=%lg, ngrid=%d, npart=%d\n",s.boxsize,s.ngrid,s.npart);
      printf("drift: [%lg  %lg  %lg]\n", s.drift[0], s.drift[1], s.drift[2]);
      printf("initial drift: [%lg  %lg  %lg]\n",s.drift_init[0],s.drift_init[1],s.drift_init[2]);


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
      pid_fname=iniparser_getstring(dict,"Rect:particle_id_file_name", "x.dat");

      fname_pinit=iniparser_getstring(dict,"Rect:init_particle_file_name", "x.dat");
      fname_pid_init=iniparser_getstring(dict,"Rect:init_particle_id_file_name","x.dat");


      raw_disp_fname=iniparser_getstring(dict,"Rect:raw_disp_field_fname", "x.dat");
      stat_disp_fname=iniparser_getstring(dict,"Rect:stat_disp_field_fname", "x.dat");

      phi_mlik_fname=iniparser_getstring(dict,"Rect:stat_phi_mlik_fname", "x.dat");
      mlik_rec_out_fname=iniparser_getstring(dict,"Rect:phi_mlik_rec_out_fname", "x.dat");

      /*-----------------------------------------------------------------------*/
      // ->> read controller <<- // 
      rc.do_rect=iniparser_getboolean(dict, "Rect:do_reconstruction", INIFALSE);
      rc.displacement_intp=iniparser_getboolean(dict, "Rect:displacement_interpolation", INIFALSE);
      rc.do_disp_perturb=iniparser_getboolean(dict, "Rect:perturbe_displacement", INIFALSE);
      //->> reconstruction type, displacement type etc. <<- //
      rc.rec_type=iniparser_getstring(dict, "Rect:reconstruction_type", "displaced");
      rc.displacement_type=iniparser_getstring(dict,"Rect:displacement_type", "backward");
      rc.displacement_order=iniparser_getstring(dict,"Rect:displacement_order", "1LPT");

      rc.displacement_tf_fname=iniparser_getstring(dict,"Rect:disp_transfunc_fname", "None");

      // ->> other controller <<- //
      do_density=iniparser_getboolean(dict, "Rect:do_density", INIFALSE);
      import_density=iniparser_getboolean(dict, "Rect:import_density", INIFALSE);

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
      // ->> testing initialization <<- //
      int do_fftw_testing;
      char *fftw_test_fname, *fftw_return_type;
      do_fftw_testing=iniparser_getboolean(dict, "Rect:do_fftw_testing", INIFALSE);
      fftw_test_fname=iniparser_getstring(dict,"Rect:fftw_test_fname", "y.dat");
      fftw_return_type=iniparser_getstring(dict,"Rect:fftw_test_return_type", "gradient");

      s.test_fname=iniparser_getstring(dict,"Rect:other_test_fname", "y.dat");

      // ->> testing for displacement/likelihood <<- //
      int do_likelihood_testing;
      char *likeli_test_fname;

      do_likelihood_testing=iniparser_getboolean(dict, "Rect:do_likelihood_testing", INIFALSE);
      likeli_test_fname=iniparser_getstring(dict,"Rect:likelihood_test_fname", "y.dat");

    /*-----------------------------------------------------------------------*/
    /*-----     End of initialization.    ------*/


    /*-----------------------------------------------------
            // ->>     Importing data     <<- //
    -----------------------------------------------------*/
    //  ->> loading particle data <<- //
    Pdata_pos *p=(Pdata_pos *)malloc(s.npart*sizeof(Pdata_pos));
    if(strcmp(rc.displacement_type, "backward_displacement")==0 ) {
      load_cita_simulation_position(particle_fname, p, s.npart);
      }
    else if(strcmp(rc.displacement_type, "likelihood_reconstruction")==0) {
      load_cita_simulation_position_pid(particle_fname, pid_fname, p, s.npart);
      }
    else{
      myerr("Simulation loading error.", FALSE);
      }
    

    // ->> initialize the boundary of particles <<- //
    s.pmin=(double *)malloc(3*sizeof(double));
    s.pmax=(double *)malloc(3*sizeof(double));
    s.dpart=(double *)malloc(3*sizeof(double));

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

      // ->> CIC density estimation, get particle boundary as well <<- //
      dmean=cic_density(p, d, s.boxsize, s.particle_mass, s.npart, s.ngrid_xyz, &s); 

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
      if(import_density==TRUE){ 
        // ->> otherwise, import density field <<- //
        main_dtype="float";
        printf("Importing the original density map directly, dtype=%s\n", main_dtype);
        load_scalar_map(oden_fname, d, s.ngrid, main_dtype);

        // ->> get particle boundary <<- //
        get_particle_boundary(p, s.boxsize, s.npart, s.ngrid_xyz, s.pmin, s.pmax, s.dpart);
        }
      else{
        printf("\n->> CAUTION: NO density map were estimated nor imported <<-\n\n"); 
	fflush(stdout);
        }
      }


    // ->> Testing <<- //
    if(do_likelihood_testing==TRUE){
      printf("Do Likelihood testing.\n"); fflush(stdout);
      //test_displacement(&s, p, d, fname_pinit, likeli_test_fname);
      //test_disp_vel_comp(&s, p, d, fname_pinit, likeli_test_fname);
      
      test_disp_direct_cal(&s, p, d, fname_pinit, likeli_test_fname);

      // ->> free <<- //
      iniparser_freedict(dict);
      myinterp_free(power); free(power);
      free(p); free(d);
      free(s.pmin); free(s.pmax); free(s.dpart);

      goto stop;
      }

    if(do_fftw_testing==TRUE) {
      printf("Do FFTW testing.\n"); fflush(stdout);
      fftw_tester(&s, d, fftw_return_type, fftw_test_fname);

      // ->> free <<- //
      iniparser_freedict(dict);
      myinterp_free(power); free(power);
      free(p); free(d);
      free(s.pmin); free(s.pmax); free(s.dpart);

      goto stop;
      }

    /*-----------------------------------------------------
         // ->>   performing Backward reconstruction   <<- //
    -----------------------------------------------------*/
    if(rc.do_rect!=TRUE) {
      printf("Do NOT perform reconstruction.\n");
      fflush(stdout); 
      goto stop; 
      }

    // ->> Obtain displacement field <<- //
    size_t ng_trim, ntrim, npart_trim;
    double bsize_trim;
    ntrim=5; //ntrim=0; 

    ng_trim=s.ngrid-2*ntrim;
    npart_trim=ng_trim*ng_trim*ng_trim;
    bsize_trim=s.boxsize*(double)ng_trim/(double)s.ngrid;

    // ->> store trim info in 'SimInfo' as well <<- //
    s.ng_trim=ng_trim; s.ntrim=ntrim;  s.npart_trim=npart_trim;
    s.bsize_trim=bsize_trim;

    //->> store trim info in 'Pdata_pos'
    cp_pdata_info(&s, p);

    float *drec, *d_disp, *d_shift;
    float *disp, *disp_lpt, *disp_mc;   //->>displacement field
    float *disp_lpt_trim, *disp_trim, *disp_mc_trim;  // ->> trimmed displacement field
    float *div, *phi, *disp_phi, *div_lpt, *phi_lpt, *disp_phi_lpt;

    //->> histogram <<- //
    //Histopar2d hist;
     

    if(rc.do_rect==TRUE){
      // ->> if do reconstruction <<- //
    
      if(strcmp(rc.displacement_type, "backward_displacement")==0 ) {
        drec=(float *)malloc(sizeof(float)*s.ngrid*s.ngrid*s.ngrid);
        d_disp=(float *)malloc(sizeof(float)*s.ngrid*s.ngrid*s.ngrid);
        d_shift=(float *)malloc(sizeof(float)*s.ngrid*s.ngrid*s.ngrid);
        
        reconstruction_partmover(&rc, &s, p, d, drec, d_disp, d_shift);

        // ->> write files <<- //
        fp=fopen(rc.rec_fname, "wb");

        fwrite(drec, sizeof(float), s.ngrid*s.ngrid*s.ngrid, fp);
        fwrite(d_disp, sizeof(float), s.ngrid*s.ngrid*s.ngrid, fp);
        fwrite(d_shift, sizeof(float), s.ngrid*s.ngrid*s.ngrid, fp);

        // ->> smooth d_disp with `inverse_gaussian' 1.-W(k,R) <<- //
	if(FALSE!=FALSE) {
          smooth_field(d_disp, s.boxsize, s.ngrid, _INVERSE_GAUSSIAN_SMOOTH_, s.smooth_R, NULL);
          fwrite(d_disp, sizeof(float), s.ngrid*s.ngrid*s.ngrid, fp);
	  }
        fclose(fp);

        // ->> free <<- //
        free(d_shift); free(d_disp); free(drec);
        }

      else if(strcmp(rc.displacement_type, "likelihood_reconstruction")==0) {
        // ->> initialization of transfer function <<- //

	// ->> allocate memory for LPT and full displacement field <<- //
        disp_lpt=(float *)fftwf_malloc(sizeof(float)*s.ngrid*s.ngrid*s.ngrid*3);
        disp=(float *)fftwf_malloc(sizeof(float)*s.ngrid*s.ngrid*s.ngrid*3);
	
        Interpar *tf=(Interpar *)malloc(3*sizeof(Interpar));
	if (transfer_func_init(tf, rc.displacement_tf_fname)!=TRUE){
          // ->> if there's no transfer function file, generate necessary data for <<- //
	  printf("Transfer function initialization failed, output displacement instead.\n"); 
	  fflush(stdout);
          load_displacement(&cp, &s, p, disp, disp_lpt, fname_pinit, fname_pid_init);
          output_real_disp_field(disp, disp_lpt, s.ngrid, raw_disp_fname);

	  fftwf_free(disp);  fftwf_free(disp_lpt);
          goto stop;
	  }

        // ->> phi max-likelihood interpolator <<- //
        Interpar *phi_mlik=(Interpar *)malloc(sizeof(Interpar));
        if(phi_mlik_init(phi_mlik, phi_mlik_fname)!=TRUE) {
	  printf("max likelihood function initialization failed, output statistic displacement field instead.\n"); 
	  fflush(stdout);

	  // ->> allocate memory <<- //
          disp_lpt_trim=(float *)fftwf_malloc(sizeof(float)*npart_trim*3);
          disp_trim=(float *)fftwf_malloc(sizeof(float)*npart_trim*3);
          disp_mc_trim=(float *)fftwf_malloc(sizeof(float)*npart_trim*3);
	  // ->> 
	  div=(float *)fftwf_malloc(sizeof(float)*npart_trim);
	  phi=(float *)fftwf_malloc(sizeof(float)*npart_trim);
	  div_lpt=(float *)fftwf_malloc(sizeof(float)*npart_trim);
	  phi_lpt=(float *)fftwf_malloc(sizeof(float)*npart_trim);
	  disp_phi=(float *)fftwf_malloc(sizeof(float)*npart_trim*3);
	  disp_phi_lpt=(float *)fftwf_malloc(sizeof(float)*npart_trim*3);
           
          load_displacement(&cp, &s, p, disp, disp_lpt, fname_pinit, fname_pid_init);

	  // ->> trim the boundary of original data <<- //
	  for(i=0; i<3; i++) {
            cubic_trim(&disp[i*s.ngrid*s.ngrid*s.ngrid], 
	               &disp_trim[i*npart_trim], s.ngrid, ntrim);
            cubic_trim(&disp_lpt[i*s.ngrid*s.ngrid*s.ngrid], 
	               &disp_lpt_trim[i*npart_trim], s.ngrid, ntrim);
	    }

          // ->> statistical separate LPT & mode-coupling term <<- //
          disp_stat_separation(&cp, disp_trim, disp_lpt_trim, disp_mc_trim, tf, 
	                       bsize_trim, npart_trim, ng_trim);

          // ->> get potential field <<- //
          potential_curlfree_vec(disp_trim, div, phi, disp_phi, bsize_trim, ng_trim);
          potential_curlfree_vec(disp_lpt_trim,div_lpt,phi_lpt,disp_phi_lpt,bsize_trim,ng_trim);
          //potential_curlfree_vec(disp_lpt_trim,div_lpt,phi_lpt,NULL,bsize_trim, ng_trim);

          // ->> output displacement field <<- //
          output_stat_disp_potential_model(disp_trim, disp_lpt_trim, disp_mc_trim, div, 
                                       phi, disp_phi, div_lpt, phi_lpt, disp_phi_lpt, 
				       s.ngrid, ng_trim, stat_disp_fname);

	  fftwf_free(disp);  fftwf_free(disp_lpt);  fftwf_free(disp_mc_trim);
	  fftwf_free(disp_trim);   fftwf_free(disp_lpt_trim);

          fftwf_free(div);  fftwf_free(phi);  fftwf_free(disp_phi); 
	  fftwf_free(div_lpt); fftwf_free(phi_lpt);

	  goto local_free;
	  }

        //->> start to potential max-likelihood reconstruction <<- //
	//->> check if I have imported position <<- //
        phi_mlik_displacement(&s, p, phi_mlik, disp, disp_phi, disp_lpt, phi, phi_lpt, 
                          ng_trim, bsize_trim, stat_disp_fname, mlik_rec_out_fname, TRUE);
          



        // ->> free all <<- //
        local_free:
        transfer_func_finalize(tf);

        }

      }




    /*-----------------------------------------------------
                       free all
    -----------------------------------------------------*/
    iniparser_freedict(dict);
    myinterp_free(power); free(power);

    free(p); free(d);
    free(s.pmin); free(s.pmax); free(s.dpart);

    //if(rc.do_rect==TRUE){
    //  free(d_shift); free(d_disp); free(drec);
    //  }


    stop:
    #ifdef _MPI_
      MPI_Finalize();
      //MPI_Barrier(MPI_COMM_WORLD);

      if(mpi_rank==0)  {
        t1=clock();
        printf("->> running time: %dsec\n\n", (int)((t1-t0)/CLOCKS_PER_SEC));
        fflush(stdout);
	}

    #else 
      t1=clock();
      printf("->> running time: %dsec\n\n", (int)((t1-t0)/CLOCKS_PER_SEC));
      fflush(stdout);
    #endif


    }






// ->> now get histogram data <<- //
/*
for(i=0; i<2; i++) {
  hist.grid[i]=100;
  hist.boundary[0][i]=200*pow(-1.,i);
  hist.boundary[1][i]=10*pow(-1.,i); }
histogram2d_init(&hist);

//histogram2d_free(&hist);
*/


