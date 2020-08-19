#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: bigden.o biglag.o calfun.o newuoa.o \
      newuob.o trsapp.o update.o main.o
	$(FC) $(FFLAGS) -o main bigden.o biglag.o calfun.o newuoa.o \
                        newuob.o trsapp.o update.o main.o
test: main
	main > output ; diff output output.txt ; rm -f output
bigden.o: bigden.f
	$(FC) $(FFLAGS) -c bigden.f
biglag.o: biglag.f
	$(FC) $(FFLAGS) -c biglag.f
calfun.o: calfun.f
	$(FC) $(FFLAGS) -c calfun.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
trsapp.o: trsapp.f
	$(FC) $(FFLAGS) -c trsapp.f
newuoa.o: newuoa.f
	$(FC) $(FFLAGS) -c newuoa.f
newuob.o: newuob.f
	$(FC) $(FFLAGS) -c newuob.f
update.o: update.f
	$(FC) $(FFLAGS) -c update.f
clean:
	rm -f *.o main
