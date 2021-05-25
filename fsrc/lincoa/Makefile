#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: calfun.o getact.o lincoa.o lincob.o \
      prelim.o qmstep.o trstep.o update.o main.o
	$(FC) $(FFLAGS) -o main calfun.o getact.o lincoa.o lincob.o \
                                prelim.o qmstep.o trstep.o update.o main.o
test: main
	./main ; rm -f *.o main
calfun.o: calfun.f90
	$(FC) $(FFLAGS) -c calfun.f90
getact.o: getact.f90
	$(FC) $(FFLAGS) -c getact.f90
lincoa.o: lincoa.f90
	$(FC) $(FFLAGS) -c lincoa.f90
lincob.o: lincob.f90
	$(FC) $(FFLAGS) -c lincob.f90
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90
prelim.o: prelim.f90
	$(FC) $(FFLAGS) -c prelim.f90
qmstep.o: qmstep.f90
	$(FC) $(FFLAGS) -c qmstep.f90
trstep.o: trstep.f90
	$(FC) $(FFLAGS) -c trstep.f90
update.o: update.f90
	$(FC) $(FFLAGS) -c update.f90
clean:
	rm -f *.o main
