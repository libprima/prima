#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: pdfomod.o bigden.o biglag.o calfun.o newuoa.o \
      newuob.o trsapp.o update.o main.o
	$(FC) $(FFLAGS) -o main pdfomod.o bigden.o biglag.o calfun.o newuoa.o \
                        newuob.o trsapp.o update.o main.o
test: main
	main > output ; diff output output.txt ; rm -f output
pdfomod.o: ../pdfomod.F
	$(FC) $(FFLAGS) -c ../pdfomod.F
bigden.o: bigden.f90
	$(FC) $(FFLAGS) -c bigden.f90
biglag.o: biglag.f90
	$(FC) $(FFLAGS) -c biglag.f90
calfun.o: calfun.f
	$(FC) $(FFLAGS) -c calfun.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
trsapp.o: trsapp.f90
	$(FC) $(FFLAGS) -c trsapp.f90
newuoa.o: newuoa.f90
	$(FC) $(FFLAGS) -c newuoa.f90
newuob.o: newuob.f
	$(FC) $(FFLAGS) -c newuob.f
update.o: update.f90
	$(FC) $(FFLAGS) -c update.f90
clean:
	rm -f *.o main
