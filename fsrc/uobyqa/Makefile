#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: lagmax.o trstep.o uobyqb.o uobyqa.o calfun.o main.o
	$(FC) $(FFLAGS) -o main lagmax.o trstep.o uobyqb.o uobyqa.o \
                                calfun.o main.o
test: main
	./main ; rm -f *.o main
calfun.o: calfun.f90
	$(FC) $(FFLAGS)  -c calfun.f90
lagmax.o: lagmax.f90
	$(FC) $(FFLAGS)  -c lagmax.f90
main.o: main.f90
	$(FC) $(FFLAGS)  -c main.f90
trstep.o: trstep.f90
	$(FC) $(FFLAGS)  -c trstep.f90
uobyqa.o: uobyqa.f90
	$(FC) $(FFLAGS)  -c uobyqa.f90
uobyqb.o: uobyqb.f90
	$(FC) $(FFLAGS)  -c uobyqb.f90
clean:
	rm -f *.o main
