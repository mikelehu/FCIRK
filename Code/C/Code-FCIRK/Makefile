#
# make : 2018-07-05
#
# Description: Makefile for compiling FCIRK project.
# System: Ubuntu 18.04 lts (64 bits)
#
#    $ make all
#
# Generate four libraries:
#     libFCIRKDOUBLE2LONG: FCIRK-DL double-long mixed-precision      
#     libFCIRKDOUBLE: FCIRK-DQ double-quad mixed-precision
#     libFCIRKLONG:   FCIRK-LQ long-quad mixed-precision
#     libFCIRKQUAD:   FCIRK-Q quadruple version.
#


# Compiler.
SYS = Linux
#CC = /usr/bin/gcc 
CC = /usr/bin/gcc 

# FLAGS.
MPFRPATH=/usr/local/include/
FLAGS=-I$(MPFRPATH) -O2 -Wall -std=c99 -fno-common -mfma

#LIBFLAGS.
LIBPATH1=/usr/local/lib/
LIBFLAGS=-L$(LIBPATH1) 

#OPENMP
OMP_FLAGS=-fopenmp

# Objects
DOUBLE2LONGOBJECTS = Common-FCIRK.o Problems.o GaussCoefficients.o Kepler.o Julia_double_FCIRK.o
DOUBLEOBJECTS = Common-FCIRK.o Problems.o GaussCoefficients.o Kepler.o Julia_double_FCIRK.o
LONGOBJECTS = Common-FCIRK.o Problems.o GaussCoefficients.o Kepler.o Julia_long_FCIRK.o
QUADOBJECTS = Common-FCIRK.o Problems.o GaussCoefficients.o Kepler.o Julia_quad_FCIRK.o

.PHONY: clean help

.SILENT: 

all :
	make clean;
	make double2long-FCIRK;
	make double-FCIRK;
	make long-FCIRK;
	make quad-FCIRK;

#
# PARALLEL: Active opem-mp parallel execution
# IOUT :    Active result in a file 
#-D DEBUG:  Active general debug
# NODEBUG
%.o: %.c def.h
	${CC} -D DEBUG -D IOUT -D PREC=$(MYPREC) -D PARALLEL -D MPFR_WANT_FLOAT128 -fPIC ${FLAGS} ${OMP_FLAGS} -c -I./  $<


# PREC=1: double-long;
double2long-FCIRK : MYPREC := 1
double2long-FCIRK : $(DOUBLE2LONGOBJECTS)
	${CC}  -shared -o libFCIRKDOUBLE2LONG.so -fPIC $(DOUBLE2LONGOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lm -lrt -lstdc++ -luuid -ldl -lquadmath -lmpfr -lgmp; \
	rm -rf *.o *~ core *tm.* *.exe 


# PREC=2: double-quad;
double-FCIRK : MYPREC := 2
double-FCIRK : $(DOUBLEOBJECTS)
	${CC}  -shared -o libFCIRKDOUBLE.so -fPIC $(DOUBLEOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lm -lrt -lstdc++ -luuid -ldl -lquadmath -lmpfr -lgmp; \
	rm -rf *.o *~ core *tm.* *.exe 

# PREC=3; long-quad;
long-FCIRK :  MYPREC := 3	
long-FCIRK :  $(LONGOBJECTS)	
	${CC}  -shared -o libFCIRKLONG.so -fPIC $(LONGOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lm -lrt -lstdc++ -luuid -ldl -lquadmath -lmpfr -lgmp; \
	rm -rf *.o *~ core *tm.* *.exe 

# PREC=11: quad-quad
quad-FCIRK :  MYPREC := 11
quad-FCIRK :  $(QUADOBJECTS)	
	${CC}  -shared -o libFCIRKQUAD.so -fPIC $(QUADOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lm -lrt -lstdc++ -luuid -ldl -lquadmath -lmpfr -lgmp;\
	rm -rf *.o *~ core *tm.* *.exe 

clean :
	rm -rf *.o *~ core *tm.* *.exe

cleanLIB :
	rm *.so

help:
	@echo "Valid targets;"
	@echo " double-FCIRK"
	@echo " long-FCIRK"
	@echo " quad-FCIRK"

