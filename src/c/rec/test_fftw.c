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

    float *d, *phi1, *phi2;
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
        // ->>           Testing FFTW           <<- //
    -----------------------------------------------------*/

    // ->> allocate memory <<- //
    phi1=(float *)fftwf_malloc(sizeof(float)*ngrid*ngrid*ngrid);
    phi2=(float *)fftwf_malloc(sizeof(float)*ngrid*ngrid*ngrid*2);

    //->> 
    int dsize, dksize, l, m, n;
    float kmin, kx, ky, kz, sin2x, sin2y, sin2z; 
    kmin=2.*pi/boxsize;

    dsize=ngrid*ngrid*ngrid*sizeof(float);
    dksize=ngrid*ngrid*(ngrid/2+1)*sizeof(fftwf_complex);

    fftwf_complex *dk1, *dk2;

    dk1=(fftwf_complex *)fftwf_malloc(dksize);
    dk2=(fftwf_complex *)fftwf_malloc(2*dksize);


    fftwf_plan pforward, pbackward_1, pbackward_2;

    pforward=fftwf_plan_dft_r2c_3d(ngrid, ngrid, ngrid, d, dk1, FFTW_ESTIMATE);
    fftwf_execute(pforward);  

    for(l=0; l<ngrid; l++)
      for(m=0; m<ngrid; m++)
        for(n=0; n<(ngrid/2+1); n++) {

          if(l<ngrid/2) kx=l*kmin;
	  else kx=(l-ngrid)*kmin;

          if(m<ngrid/2) ky=m*kmin;
	  else ky=(m-ngrid)*kmin;

          if(n<ngrid/2) kz=n*kmin;
	  else kz=(n-ngrid)*kmin;

          sin2x = 4.*sin(kx/2.)*sin(kx/2.);
          sin2y = 4.*sin(ky/2.)*sin(ky/2.);
          sin2z = 4.*sin(kz/2.)*sin(kz/2.);


          if ((l==0) && (m==0) && (n==0)) greens = 0.;
          else greens = -1./(sin2x+sin2y+sin2z);

          // ->>  assign dk1, dk2 <<- //
          ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]*=greens;
          ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]*=greens;


          ArrayAccess3D_n3(dk2, ngrid, ngrid, (ngrid/2+1), l, m, n)[0]=
	    ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[0];

          ArrayAccess3D_n3(dk2, ngrid, ngrid, (ngrid/2+1), l, m, n)[1]=
	    ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[1];

          ArrayAccess4D_n4(dk2, 2, ngrid, ngrid, (ngrid/2+1), 1, l, m, n)[0]=
	    ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[0];

          ArrayAccess4D_n4(dk2, 2, ngrid, ngrid, (ngrid/2+1), 1, l, m, n)[1]=
	    ArrayAccess3D_n3(dk1, ngrid, ngrid, (ngrid/2+1), l, m, n)[1];

	  }


    /* ->> find the inverse FFT of phi <<- */
    pbackward_1 = fftwf_plan_dft_c2r_3d(ngrid, ngrid, ngrid, dk1, phi1, FFTW_ESTIMATE);
    fftwf_execute(pbackward_1);

    //->> advanced FFTW <<- //

    int rank, howmany, *ndim, idist, odist, istride, ostride, *inembed, *onembed;
    ndim=(int *)malloc(3*sizeof(int));
    ndim[0]=ngrid; ndim[1]=ngrid; ndim[2]=ngrid;

    rank=3;
    howmany=2;
    idist=ngrid*ngrid*(ngrid/2+1);
    odist=ngrid*ngrid*ngrid;
    istride=1; ostride=1;
    inembed=ndim; onembed=ndim;

    pbackward_2=fftwf_plan_many_dft_c2r(rank, ndim, howmany, dk2, inembed, istride, idist, phi2, onembed, ostride, odist, FFTW_ESTIMATE);
    
    fftwf_execute(pbackward_2);



    // ->> destroy fftw_plan <<- //
    fftwf_destroy_plan(pforward);
    fftwf_destroy_plan(pbackward_1);
    fftwf_destroy_plan(pbackward_2);

    printf("->> FFTW Done.\n");




    int _write_testfile_=TRUE;
    if(_write_testfile_){
      fp=fopen(test_fname, "wb");
      
      fwrite(phi1, sizeof(float), ngrid*ngrid*ngrid, fp);
      fwrite(phi2, sizeof(float), ngrid*ngrid*ngrid*2, fp);

      fclose(fp);
      }


    /*-----------------------------------------------------
                       free all
    -----------------------------------------------------*/
    iniparser_freedict(dict);
    free(p); free(d);

    fftwf_free(phi1);
    fftwf_free(phi2);


    }


