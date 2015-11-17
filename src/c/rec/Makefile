#---------------------------#
DO_MPI= no
#DO_MPI= yes


ifeq ($(DO_MPI), yes)
  CC= mpicc
else
  CC= gcc
endif
MPICC= mpicc
FC= ifort
#---------------------------#



# C main() call fortran subroutine, for ifort.
COMFFLAGS = -nofor_main    
FLIB = -L/usr/local/lib  -lm   

#-----------------------------------------------------------------#
#------------------------- My NOTEBOOKS---------------------------#
#-----------------------------------------------------------------#

ifeq ($(MACHINE),mylaptop)
    MYLIBPATH= /Users/harrygun/workspace/lib/mylib
    INIPATH = /Users/harrygun/workspace/lib/iniparser
    MPIPATH = /Users/harrygun/software/work/mpich2/mpich2-1.2.1p1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/harrygun/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -arch i386 -arch x86_64 -I /Users/harrygun/workspace/lib/gsl/include \
		 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)

	#CFLAGS = -g -I /Users/harrygun/workspace/lib/gsl/include \
	#	 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)
    endif

endif

#---- gwln ----
ifeq ($(MACHINE),gwln)
    MYLIBPATH= /home/wangxin/workspace/lib/mylib
    INIPATH = /home/wangxin/workspace/lib/iniparser/$(HOSTNAME)
    MPIPATH = /home/wangxin/software/mpich2/gwln1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)/$(HOSTNAME)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D_MPI_ -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -fPIC -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif

#---- HHPC ----
ifeq ($(MACHINE),HHPC)
    MYLIBPATH= /home/wangxin/lib/mylib
    INIPATH = /home/wangxin/lib/iniparser
    MPIPATH = /usr/mpi/intel/openmpi-1.4.3

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D_MPI_ -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif


#---- CITA VOLE ----#
ifeq ($(MACHINE),cita)
    MYLIBPATH= /home/xwang/workspace/lib/mylib
    INIPATH = /home/xwang/workspace/lib/iniparser
    MPIPATH = 

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

endif









#################################################################
#MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


LIB  = $(INILIB) $(GSLIB) $(MYLIB)
#lib  = $(inilib) $(gslib) 




.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<





# workspace path
WORKSPACE = ../../workspace
AR=ar
ARFLAGS=-r
RANLIB=ranlib


default: biasrpt
#default: rsdrpt_cont
#default: rsdrpt rsdrpt_cont
#default: rsdmatt


#--------------------------------------------------#
rn_integration: rn_integration.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)
#------------ General Bias ------------#

libkernel.a: kernel.o
	@($(AR) $(ARFLAGS) $@ $^)
	@($(RANLIB) $@)

rsdmatt: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o biasrpt_gal_rsd.o driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)


biasrpt: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o biasrpt_gal_rsd.o driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

rsdrpt: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o biasrpt_gal_rsd.o driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

rsdrpt_cont: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o biasrpt_gal_rsd.o driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)




biasrpt_ge1: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o  driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)
gbiasrpt: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o  driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)
GE2_biasrpt: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o biasrpt_test.o \
         biasrpt_b.o utilities.o biasrpt_G1.o biasrpt_G2.o fnl_utils.o mpinit.o\
	 fnl_rpt.o biasrpt_galax.o rsd_Gn.o  driver.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

#--------------------------------------------------#
fnl_nltrans: fnl_utils.o kernel.o fnl_spt.o fnl_rpt.o utilities.o driver_fnl.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)
#--------------------------------------------------#
#--------------------------------------------------#
#--------------------------------------------------#
#--------------------------------------------------#
logept_resum: myradmul_tmp.o integration_tmp.o  propagator.o kernel.o \
             utilities.o ept_log_prop2p.o ept_log_resum.o driver_ept_log_resum.o 
	$(CC) $^ -o $(WORKSPACE)/$@  $(LIB)
#--------------------------------------------------#



logbias:  utilities.o log_bias.o  ept_driver.o 
	$(CC) $^  -o $(WORKSPACE)/$@  $(LIB)

pdensity: propagator.o kernel.o utilities.o pdensity.o pdensity_driver.o
	$(CC) $^ -o $(WORKSPACE)/$@  $(LIB)


logept:  ept_log.o ept_log_2lp.o ept_driver.o 
	$(CC) $^  -o $(WORKSPACE)/$@  $(LIB)


loglpt:   lpt_log.o lpt_driver.o test.o
	$(CC) $^  -o $(WORKSPACE)/$@  $(LIB)


test:   lpt_log.o test_driver.o test.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

test_time: myradmul_tmp.o integration_tmp.o utilities.o propagator.o kernel.o \
           ept_log_resum.o ept_log_prop2p.o timer_test.o 
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

omp_test: omp_test.o
	$(CC)  $^ -fopenmp -o $(WORKSPACE)/$@  $(LIB)
#--------------------------------------------------

mpi_test: mpi_test.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)


#--------------------------------------------------
#            primordial nongaussianity           # 
#--------------------------------------------------
fnl_tmp:  fnl_tmp.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)

#--------------------------------------------------
clean:
	rm -f *.o *.c~



#-----------------------------------------
tmp: infsum.o tmp.o
	$(CC)   $^  -o $(WORKSPACE)/$@  $(LIB)

