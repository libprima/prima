################ Makefile by Zaikun Zhang (www.zhangzk.net) ###############

#  please change the variables FC and FFLAGS to match your fortran compiler

# FC = g95
# FC = f77
# FC = ifort
FC = gfortran
FFLAGS = -g

test: cobyla.o cobylb.o trstlp.o calcfc.o test.o
	$(FC) $(FFLAGS) -o test cobyla.o cobylb.o trstlp.o calcfc.o test.o
test: test
	./test ; rm -f *.o test
cobyla.o: cobyla.f
	$(FC) $(FFLAGS) -c cobyla.f
cobylb.o: cobylb.f
	$(FC) $(FFLAGS) -c cobylb.f
trstlp.o: trstlp.f
	$(FC) $(FFLAGS) -c trstlp.f
calcfc.o: calcfc.f90
	$(FC) $(FFLAGS) -c calcfc.f90
test.o: test.f
	$(FC) $(FFLAGS) -c test.f90
clean:
	rm -f *.o test
