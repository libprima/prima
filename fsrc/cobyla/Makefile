################ Makefile by Zaikun Zhang (www.zhangzk.net) ###############

#  please change the variables FC and FFLAGS to match your fortran compiler

# FC = g95
# FC = f77
# FC = ifort
FC = gfortran
FFLAGS = -g

main: cobyla.o cobylb.o trstlp.o calcfc.o main.o
	$(FC) $(FFLAGS) -o main cobyla.o cobylb.o trstlp.o calcfc.o main.o
test: main
	./main ; rm -f *.o main
cobyla.o: cobyla.f90
	$(FC) $(FFLAGS) -c cobyla.f90
cobylb.o: cobylb.f90
	$(FC) $(FFLAGS) -c cobylb.f90
trstlp.o: trstlp.f90
	$(FC) $(FFLAGS) -c trstlp.f90
calcfc.o: calcfc.f90
	$(FC) $(FFLAGS) -c calcfc.f90
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90
clean:
	rm -f *.o main
