buildir = .
include $(buildir)/make.config

objects = dvr input.o T.o key.o points.o V.o sort.o fprint.o crd.o fname.o gmat.o

all : $(objects)
	cp dvr $(bindir)/dvr

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

key.o : key.f90
	$(FC) $(FFLAGS) -c key.f90

fname.o : fname.f90
	$(FC) $(FFLAGS) -c fname.f90

crd.o : crd.f90
	$(FC) $(FFLAGS) -c crd.f90

sort.o : sort.f90
	$(FC) $(FFLAGS) -c sort.f90

gmat.o : gmat.f90
	$(FC) $(FFLAGS) -c gmat.f90

linal.o : linal.f90 sort.o
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 sort.o

points.o : points.f90 key.o
	$(FC) $(FFLAGS) -c points.f90 key.o

fprint.o : fprint.f90 key.o fname.o
	$(FC) $(FFLAGS) -c fprint.f90 key.o fname.o

V.o : V.f90 key.o
	$(FC) $(FFLAGS) -c V.f90 key.o

T.o : T.f90 key.o 
	$(FC) $(FFLAGS) -c T.f90 key.o

dvr : dvr.f90 input.o T.o key.o points.o V.o linal.o sort.o fprint.o crd.o fname.o gmat.o
	$(FC) $(FFLAGS) $(LINAL) -o dvr dvr.f90 input.o T.o key.o points.o V.o linal.o sort.o fprint.o crd.o fname.o gmat.o
