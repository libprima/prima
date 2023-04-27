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
	main > output ; diff output output.txt ; rm -f output
cobyla.o: cobyla.f
	$(FC) $(FFLAGS) -c cobyla.f
cobylb.o: cobylb.f
	$(FC) $(FFLAGS) -c cobylb.f
trstlp.o: trstlp.f
	$(FC) $(FFLAGS) -c trstlp.f
calcfc.o: calcfc.f
	$(FC) $(FFLAGS) -c calcfc.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
clean:
	rm -f *.o main
