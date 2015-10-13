ifeq ($(MACHINE),mylaptop)
	cflag= "-arch i386 -arch x86_64 -I /Users/harrygun/workspace/software/works/cuba/include    -I/Users/harrygun/workspace/lib/gsl/include   -I ../c/"
	libflag= "-L /Users/harrygun/workspace/software/works/cuba/lib -L /Users/harrygun/workspace/lib/mylib/script/clib  -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas  -L ./ -L ../c/"
endif

ifeq ($(MACHINE),gwln)
	cflag="-I /home/wangxin/workspace/software/cuba/include  -I/home/wangxin/workspace/lib/include  -I ../c/"
	libflag="-L /home/wangxin/workspace/software/cuba/lib -L ./ -L/home/wangxin/workspace/lib/lib -lgsl -lgslcblas -lm -L ../c/"

endif

ifeq ($(MACHINE),cita)
	cflag="-I /home/xwang/workspace/lib/gsl/include  -I ../c/"
	libflag="-L ./ -L /home/xwang/workspace/lib/gsl/lib -lgsl -lgslcblas -lm -L ../c/"

endif

all:
	CFLAGS=$(cflag)  LDFLAGS=$(libflag) \
	python setup.py build_ext --inplace 

clean:
	@echo Cleaning
	@rm -f *.c *.o *.so *~
	@rm -rf build
