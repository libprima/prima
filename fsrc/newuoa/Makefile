#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: bigden.o biglag.o calfun.o newuoa.o \
      newuob.o trsapp.o update.o main.o
	$(FC) $(FFLAGS) -o main bigden.o biglag.o calfun.o newuoa.o \
                        newuob.o trsapp.o update.o main.o
test: main
	main > output ; diff output output.txt ; rm -f output
bigden.o: bigden.f90
	$(FC) $(FFLAGS) -c bigden.f90
biglag.o: biglag.f90
	$(FC) $(FFLAGS) -c biglag.f90
calfun.o: calfun.f90
	$(FC) $(FFLAGS) -c calfun.f90
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90
trsapp.o: trsapp.f90
	$(FC) $(FFLAGS) -c trsapp.f90
newuoa.o: newuoa.f90
	$(FC) $(FFLAGS) -c newuoa.f90
newuob.o: newuob.f90
	$(FC) $(FFLAGS) -c newuob.f90
update.o: update.f90
	$(FC) $(FFLAGS) -c update.f90
clean:
	rm -f *.o main
