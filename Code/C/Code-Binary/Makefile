#
#
#    2018-07-16
#    Description: Makefile for compiling CBinFileWRA.c
#    System: Ubuntu 18.04 lts (64 bits)
#
#    $ make all
#
#    Generate library:
#       libCBinFilesWRA
#


# Compiler.
SYS = Linux
CC = /usr/bin/gcc

# FLAGS.
MPFRPATH=/usr/local/include/
FLAGS=-I$(MPFRPATH) -O2 -Wall -std=c99 -fno-common -mfma

#LIBFLAGS.
LIBPATH1=/usr/local/lib/
LIBFLAGS=-L$(LIBPATH1)

# Objects
OBJECTS = CBinFilesWRA.o

.PHONY: clean help 

.SILENT: 

all : 
	make clean;
	make fun;

#-D DEBUG:  Active general debug
%.o: %.c 
	${CC} -D NODEBUG -D MPFR_WANT_FLOAT128 -fPIC ${FLAGS} -c -I./  $<

fun: $(OBJECTS) 
	${CC}  -shared -o libCBinFilesWRA.so -fPIC $(OBJECTS) $(LIBFLAGS) -lmpfr -lgmp 

.PHONY: clean

clean :
	rm -rf *.o *.so *~ core *tm.* *.exe

help:
	@echo "Valid targets;"
	@echo " clean"
	@echo " all" 
   



