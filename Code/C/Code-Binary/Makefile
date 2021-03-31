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

#MPFRPATH=/usr/local/include/
#FLAGS=-I$(MPFRPATH) -O2 -Wall -std=c99 -fno-common -mfma

INCLUDE  = /usr/local/include 
CFLAGS   = $(addprefix -I,$(INCLUDE)) -O2 -Wall -std=c99 
CFLAGS   += -shared -fno-common -mfma -fPIC -DMPFR_WANT_FLOAT128

#LIBFLAGS.

#LIBPATH1=/usr/local/lib/
#LIBFLAGS=-L$(LIBPATH1)

LIBPATH  = /usr/local/lib
LDFLAGS  = $(addprefix -L,$(LIBPATH))
LDLIBS   = -lmpfr -lgmp

# Objects

#OBJECTS = CBinFilesWRA.o

SOURCES  = $(wildcard *.c)
OBJECTS  = $(patsubst %.c,%.o,$(SOURCES))

.PHONY: clean help 

.SILENT: 

all : 
	make clean;
	make libCBinary.so;

#-D DEBUG:  Active general debug
#%.o: %.c 
#	${CC} -D NODEBUG -D MPFR_WANT_FLOAT128 -fPIC ${FLAGS} -c -I./  $<

#fun: $(OBJECTS) 
#	${CC}  -shared -o libCBinFilesWRA.so -fPIC $(OBJECTS) $(LIBFLAGS) -lmpfr -lgmp 

libCBinary.so: $(OBJECTS) 
	$(CC) $(CFLAGS) $(LDFLAGS) -o$@ $^ $(LDLIBS)

.PHONY: clean

clean :
	rm -rf *.o *.so *~ core *tm.* *.exe

help:
	@echo "Valid targets;"
	@echo " clean"
	@echo " all" 
   


