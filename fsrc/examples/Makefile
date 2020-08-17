# This is the Makefile for the tests.
# Please change the variables FC and FFLAGS to match your fortran compiler.

FC = gfortran -Wall -Wextra -fcheck=all -pedantic
#FC = ifort -warn all -check all
#FC = nagfor -C 
#FC = flang -Wall -Wextra -Mstandard 
#FC = sunf95 -w3 -u -U
#FC = af95 -m1
#FC = pgfortran -Minform=warn
FFLAGS = -g 
COMMON = ../common
NEWUOA = ../newuoa

test_newuoa: newuoa_exmp
	./newuoa_exmp ; rm -f *.o *.mod newuoa_exmp *.dbg
 
# N.B.:
# The .F90 and .f90 files will be compiled in the enumeration order of
# the .o files. The order matters, because the compilation of each .o file 
# depends on the .o files (and the corresponding .mod files) preceding it.
newuoa_exmp: newuoa_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o lina.o output.o pintrf.o \
	preproc.o initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o 
	$(FC) $(FFLAGS) -o newuoa_exmp newuoa_exmp.f90 *.o

# Compile the Fortran code providing generic modules
consts.o: $(COMMON)/consts.F90
	$(FC) $(FFLAGS) -c $(COMMON)/consts.F90
info.o: $(COMMON)/info.F90
	$(FC) $(FFLAGS) -c $(COMMON)/info.F90
debug.o: $(COMMON)/debug.F90
	$(FC) $(FFLAGS) -c $(COMMON)/debug.F90
memory.o: $(COMMON)/memory.F90
	$(FC) $(FFLAGS) -c $(COMMON)/memory.F90
infnan.o: $(COMMON)/infnan.F90
	$(FC) $(FFLAGS) -c $(COMMON)/infnan.F90
lina.o: $(COMMON)/lina.F90
	$(FC) $(FFLAGS) -c $(COMMON)/lina.F90
output.o: $(COMMON)/output.F90
	$(FC) $(FFLAGS) -c $(COMMON)/output.F90
pintrf.o: $(COMMON)/pintrf.F90
	$(FC) $(FFLAGS) -c $(COMMON)/pintrf.F90

# Compile the Fortran code providing solver-specific modules
preproc.o: $(NEWUOA)/preproc.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/preproc.f90
initialize.o: $(NEWUOA)/initialize.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/initialize.f90
trustregion.o: $(NEWUOA)/trustregion.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/trustregion.f90
vlagbeta.o: $(NEWUOA)/vlagbeta.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/vlagbeta.f90
geometry.o: $(NEWUOA)/geometry.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/geometry.f90
shiftbase.o: $(NEWUOA)/shiftbase.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/shiftbase.f90
update.o: $(NEWUOA)/update.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/update.f90
newuob.o: $(NEWUOA)/newuob.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/newuob.f90
newuoa.o: $(NEWUOA)/newuoa.f90
	$(FC) $(FFLAGS) -c $(NEWUOA)/newuoa.f90 

clean:
	rm -f *.o *.mod newuoa_exmp *.dbg
