#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: altmov.o bobyqa.o bobyqb.o calfun.o \
      prelim.o rescue.o trsbox.o update.o main.o
	$(FC) $(FFLAGS) -o main altmov.o bobyqa.o bobyqb.o calfun.o \
                           main.o prelim.o rescue.o trsbox.o update.o
test: main
	main > output ; diff output output.txt ; rm -f output 
altmov.o: altmov.f
	$(FC) $(FFLAGS) -c altmov.f
bobyqa.o: bobyqa.f
	$(FC) $(FFLAGS) -c bobyqa.f
bobyqb.o: bobyqb.f
	$(FC) $(FFLAGS) -c bobyqb.f
calfun.o: calfun.f
	$(FC) $(FFLAGS) -c calfun.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
prelim.o: prelim.f
	$(FC) $(FFLAGS) -c prelim.f
rescue.o: rescue.f
	$(FC) $(FFLAGS) -c rescue.f
trsbox.o: trsbox.f
	$(FC) $(FFLAGS) -c trsbox.f
update.o: update.f
	$(FC) $(FFLAGS) -c update.f
clean:
	rm -f *.o main
