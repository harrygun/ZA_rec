  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <fftw3.h>

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





  int main( int argc, char *argv[])  {

      int debug= 20, i, j, k, do_MPI=99;
      char *ini_name;

      if(argc!=2)  //
        myerr("Input parameters file is needed.", FALSE);

        ini_name=argv[1];

      do_MPI = FALSE;

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
      char *smooth_type, *particle_fname, *droot, *plin_name, *oden_fname, *test_fname, *main_dtype;
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

      test_fname=iniparser_getstring(dict,"Rect:test_fname", "y.dat");

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

    float *d, *phi, *phi_ij;
    double particle_mass, rhom_, dmean;
    int ngrid_xyz[3];

    for(i=0; i<3; i++)
      ngrid_xyz[i]=ngrid;
    d=(float *)malloc(sizeof(float)*ngrid_xyz[0]*ngrid_xyz[1]*ngrid_xyz[2]);

    // ->> if do CIC density estimation <<- //
    if (do_density==TRUE) {
      rhom_=get_rhom(&cp, cp.z);
      particle_mass=part_mass(&cp, cp.z, boxsize, ngrid);
      printf("\nmean density=%lg (M_star)*(h/Mpc)^3\nparticle mass resolution=%lg\n\n", rhom_, particle_mass);

      // ->> CIC density estimation <<- //
      dmean=cic_density(p, d, boxsize, particle_mass, npart, ngrid_xyz); 

      // ->> save density field in file <<- //
      if(save_odensity==TRUE) { //
        fp=fopen(oden_fname, "wb");
        for(i=0; i<ngrid; i++)
          for(j=0; j<ngrid; j++)
            for(k=0; k<ngrid; k++)
              fwrite(&ArrayAccess3D(d, ngrid, i, j, k), sizeof(float), 1, fp);
        fclose(fp);
        }

      }
    else {
      // ->> otherwise, import density field <<- //
      main_dtype="float";
      printf("Importing the original density map directly, dtype=%s\n", main_dtype);
      load_scalar_map(oden_fname, d, ngrid, main_dtype);
      }



    /*-----------------------------------------------------
            // ->>     Testing FFTW     <<- //
    -----------------------------------------------------*/


    // ->> allocate memory <<- //
    phi=(float *)fftwf_malloc(sizeof(float)*ngrid_xyz[0]*ngrid_xyz[1]*ngrid_xyz[2]);
    phi_ij=(float *)fftwf_malloc(sizeof(float)*ngrid_xyz[0]*ngrid_xyz[1]*ngrid_xyz[2]*10);

    poisson_solver_float(d, phi, phi_i, phi_ij, boxsize, ngrid, cp.flg[0], cp.R, fft_return_type);










    printf("Done.\n");














    int _write_testfile_=TRUE;
    if(_write_testfile_){
      fp=fopen(test_fname, "wb");
      
      fwrite(phi, sizeof(float), ngrid*ngrid*ngrid, fp);
      fwrite(phi_ij, sizeof(float), ngrid*ngrid*ngrid*10, fp);

      fclose(fp);
      }


    /*-----------------------------------------------------
                       free all
    -----------------------------------------------------*/
    iniparser_freedict(dict);
    free(p); free(d);
    fftwf_free(phi);

    if(do_grad) fftwf_free(phi_i);
    if(do_hess) fftwf_free(phi_ij);

    }


