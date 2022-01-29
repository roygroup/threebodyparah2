CC=g++
CFLAGS=-O3

F90=gfortran
FFLAGS=-O3
ffile=f90RKHS

options=-std=c++11 

LDLIBS=-lgfortran 

OBJS= extrapolate.o jpoint.o callRKHS.o RKHS.o region.o trapezoid.o utils.o rkhs.o


all: threebodyparah2 clean

threebodyparah2: main.o $(OBJS)
	$(CC) -o $@ main.o $(OBJS) $(LDLIBS)

extrapolate.o:   extrapolate.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

jpoint.o:   jpoint.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

region.o:   region.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

trapezoid.o:   trapezoid.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

utils.o:   utils.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

rkhs.o:   rkhs.cpp
	$(CC) $(options) -c $(CFLAGS) -o $@ $*.cpp

callRKHS.o: $(ffile)/callRKHS.f90 RKHS.o
	$(F90) -c $(FFLAGS) -o $@ $(ffile)/callRKHS.f90

RKHS.o: $(ffile)/RKHS.f90
	$(F90) -c $(FFLAGS) -o $@ $(ffile)/RKHS.f90


.PHONY: clean
clean:
	rm *.mod *.o
