buildir = .
include $(buildir)/make.config

objects = dvr input.o T.o key.o points.o V.o sort.o wave.o crd.o

all : $(objects)
	cp dvr $(bindir)/dvr

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

key.o : key.f90
	$(FC) $(FFLAGS) -c key.f90

crd.o : crd.f90
	$(FC) $(FFLAGS) -c crd.f90

sort.o : sort.f90
	$(FC) $(FFLAGS) -c sort.f90

linal.o : linal.f90 sort.o
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 sort.o

points.o : points.f90 key.o
	$(FC) $(FFLAGS) -c points.f90 key.o

wave.o : wave.f90 key.o
	$(FC) $(FFLAGS) -c wave.f90 key.o

T.o : T.f90 key.o
	$(FC) $(FFLAGS) -c T.f90 key.o

V.o : V.f90 key.o
	$(FC) $(FFLAGS) -c V.f90 key.o

dvr : dvr.f90 input.o T.o key.o points.o V.o linal.o sort.o wave.o crd.o
	$(FC) $(FFLAGS) $(LINAL) -o dvr dvr.f90 input.o T.o key.o points.o V.o linal.o sort.o wave.o crd.o
