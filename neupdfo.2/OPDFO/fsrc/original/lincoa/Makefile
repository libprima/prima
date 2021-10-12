#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: calfun.o getact.o lincoa.o lincob.o \
      prelim.o qmstep.o trstep.o update.o main.o
	$(FC) $(FFLAGS) -o main calfun.o getact.o lincoa.o lincob.o \
                                prelim.o qmstep.o trstep.o update.o main.o 
test: main
	main > output ; diff output output.txt ; rm -f output
calfun.o: calfun.f
	$(FC) $(FFLAGS) -c calfun.f
getact.o: getact.f
	$(FC) $(FFLAGS) -c getact.f
lincoa.o: lincoa.f
	$(FC) $(FFLAGS) -c lincoa.f
lincob.o: lincob.f
	$(FC) $(FFLAGS) -c lincob.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
prelim.o: prelim.f
	$(FC) $(FFLAGS) -c prelim.f
qmstep.o: qmstep.f
	$(FC) $(FFLAGS) -c qmstep.f
trstep.o: trstep.f
	$(FC) $(FFLAGS) -c trstep.f
update.o: update.f
	$(FC) $(FFLAGS) -c update.f
clean:
	rm -f *.o main
