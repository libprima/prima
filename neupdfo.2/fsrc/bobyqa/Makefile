#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: altmov.o bobyqa.o bobyqb.o calfun.o \
      prelim.o rescue.o trsbox.o update.o main.o
	$(FC) $(FFLAGS) -o main altmov.o bobyqa.o bobyqb.o calfun.o \
                           main.o prelim.o rescue.o trsbox.o update.o
test: main
	./main ; rm -f *.o main
altmov.o: altmov.f90
	$(FC) $(FFLAGS) -c altmov.f90
bobyqa.o: bobyqa.f90
	$(FC) $(FFLAGS) -c bobyqa.f90
bobyqb.o: bobyqb.f90
	$(FC) $(FFLAGS) -c bobyqb.f90
calfun.o: calfun.f90
	$(FC) $(FFLAGS) -c calfun.f90
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90
prelim.o: prelim.f90
	$(FC) $(FFLAGS) -c prelim.f90
rescue.o: rescue.f90
	$(FC) $(FFLAGS) -c rescue.f90
trsbox.o: trsbox.f90
	$(FC) $(FFLAGS) -c trsbox.f90
update.o: update.f90
	$(FC) $(FFLAGS) -c update.f90
clean:
	rm -f *.o main
