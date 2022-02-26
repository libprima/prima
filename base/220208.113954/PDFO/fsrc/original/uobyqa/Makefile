#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: lagmax.o trstep.o uobyqb.o uobyqa.o calfun.o main.o
	$(FC) $(FFLAGS) -o main lagmax.o trstep.o uobyqb.o uobyqa.o \
                                calfun.o main.o
test: main
	main > output ; diff output output.txt ; rm -f output
calfun.o: calfun.f
	$(FC) $(FFLAGS)  -c calfun.f
lagmax.o: lagmax.f
	$(FC) $(FFLAGS)  -c lagmax.f
main.o: main.f
	$(FC) $(FFLAGS)  -c main.f
trstep.o: trstep.f
	$(FC) $(FFLAGS)  -c trstep.f
uobyqa.o: uobyqa.f
	$(FC) $(FFLAGS)  -c uobyqa.f
uobyqb.o: uobyqb.f
	$(FC) $(FFLAGS)  -c uobyqb.f
clean:
	rm -f *.o main
