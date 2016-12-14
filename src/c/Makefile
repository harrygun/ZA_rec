#---------------------------#
DO_MPI= no
#DO_MPI= yes

#-->
##DO_OMP=yes
DO_OMP=no

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
	CFLAGS = -g -arch x86_64 -D _MPI_ -I /Users/harrygun/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
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
    #MPIPATH = 
    #FFTWPATH = 

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include
    #FFTWINCLUDE = $(FFTWPATH)/

    GENLIB   = lmylib
    FISHLIB  = lfish

    #FFTWLIB = -lfftw3 -lm
    #FFTWLIB = -lfftw3f -lfftw3 -lm

    GSLIB = -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB)  $(FFTWLIB)


    ifeq ($(DO_OMP),yes)
        CFLAGS = -fopenmp -g  -D _OMP_ -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
        FFTWLIB =  -lfftw3f_threads -lfftw3_threads -lfftw3f -lfftw3 -lm

    else
        CFLAGS = -g -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
        FFTWLIB = -lfftw3f -lfftw3 -lm
    endif

endif

#----------------------#
ifeq ($(MACHINE),scinet-gpc)
    MYLIBPATH= /home/p/pen/wangxin/workspace/lib/mylib
    INIPATH = /home/p/pen/wangxin/workspace/lib/iniparser
    #MPIPATH = 
    #FFTWPATH = 

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include
    #FFTWINCLUDE = $(FFTWPATH)/

    GENLIB   = lmylib
    FISHLIB  = lfish

    #FFTWLIB = -lfftw3 -lm
    #FFTWLIB = -lfftw3f -lfftw3 -lm

    GSLIB = -L/home/p/pen/wangxin/workspace/lib/lib/ -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB)  $(FFTWLIB)


    ifeq ($(DO_OMP),yes)
        CFLAGS = -fopenmp -g  -D _OMP_ -I /home/p/pen/wangxin/workspace/lib/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
        FFTWLIB =  -lfftw3f_threads -lfftw3_threads -lfftw3f -lfftw3 -lm

    else
        CFLAGS = -g -I /home/p/pen/wangxin/workspace/lib/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
        FFTWLIB = -lfftw3f -lfftw3 -lm
    endif

endif






#################################################################
#MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


LIB  = $(INILIB)  $(MYLIB) $(GSLIB)
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


default: rect test_fftw test_field
#default: rect



#--------------------------------------------------#
#rect: driver.o io.o cic.o poisson.o partmove.o misc.o za_reconstruction.o

rect: driver.o io.o cic.o poisson.o partmove.o misc.o reconstruction_partmoving.o backward_displacement.o fourier.o stat_model.o displacement.o test_likelihood.o stat_model.o maxlike_phi_rec.o
	$(CC) $(CFLAGS)  $^ -o $(WORKSPACE)/$@  $(LIB)

test_fftw: test_fftw.o io.o cic.o
	$(CC) $(CFLAGS)  $^ -o $(WORKSPACE)/$@  $(LIB)

test_field: test_driver.o io.o cic.o poisson.o misc.o
	$(CC) $(CFLAGS)  $^ -o $(WORKSPACE)/$@  $(LIB)


#--------------------------------------------------
clean:
	rm -f *.o *.c~




