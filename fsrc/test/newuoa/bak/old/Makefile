FSTD=2003
FS=f03

GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Wuninitialized -Wmaybe-uninitialized  -Wsurprising -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow #-Warray-temporaries
IFORT=ifort -c -g -warn all -check all  -fpe0 -debug extended -fimplicit-none #-names as_is
G95=g95 -c -g -Wall -Wextra -fbounds-check -fimplicit-none
NFORT=nagfor -c -g -u -C -Warn=allocation -Warn=constant_coindexing -Warn=reallocation -Warn=subnormal -colour=error:red,warn:magenta,info:yellow
PFORT=pgfortran  -c -g -Minform=warn
AFORT=af95  -c -g -m1 -en
SFORT=sunf95  -c -g -w3 -u -ansi -U
FLANG=flang -c -Wall -Wextra -Mstandard

CURRENTDIR=$(PWD)
NEUPDFO=$(CURRENTDIR)/../../
INTERFORMDIR=../intersection_form/newuoa/
MEXDIR= $(PWD)/../../matlab/mex_gateways/
INTERSECTIONMEXDIR= $(PWD)/../../matlab/mex_gateways/intersection_form/

COMFILES=../common/consts.F90 ../common/debug.F90 ../common/memory.F90 ../common/info.F90 ../common/infnan.F90 ../common/lina.F90 ../common/output.F90 ../common/pintrf.F90
FILES= $(COMFILES) preproc.f90 initialize.f90 trustregion.f90 vlagbeta.f90 geometry.f90 update.f90 shiftbase.f90 newuob.f90 newuoa.f90

M90FILES= $(FILES) $(MEXDIR)/fmxapi.F90 $(MEXDIR)/prob.F90 $(MEXDIR)/classical/fmxcl.F90 $(MEXDIR)/gethuge.F90 $(MEXDIR)/newuoa-interface.F90 $(MEXDIR)/classical/newuoa-interface.F90 -I/home/zaikunzhang/local/MATLAB/R2020a/extern/include/ -I$(NEUPDFO)

MFILES= ../common/consts.F ../common/debug.F ../common/memory.F ../common/info.F ../common/infnan.F ../common/lina.F ../common/output.F ../common/pintrf.F preproc.f initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f $(INTERSECTIONMEXDIR)/fmxapi.F $(INTERSECTIONMEXDIR)/prob.F $(INTERSECTIONMEXDIR)/classical/fmxcl.F $(INTERSECTIONMEXDIR)/gethuge.F $(INTERSECTIONMEXDIR)/newuoa-interface.F $(INTERSECTIONMEXDIR)/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2020a/extern/include/ -I$(NEUPDFO)


verify: gverify iverify nverify pverify averify fverify sverify 9verify clean
verifym: gverifym iverifym nverifym pverifym averifym fverifym sverifym 9verifym clean

lint:
	#flint -I../common/ -a -f -g -m -p -3 -e $(FILES)
	flint -I../common/  $(FILES)

fverify:
	$(FLANG) $(FILES) -ffree-form -std=f$(FSTD)

gverify:
	$(GFORT)  $(FILES) -std=f$(FSTD) -pedantic -ffree-form

9verify:
	#$(G95)  $(FILES) -std=F -pedantic -ffree-form
	$(G95)  $(FILES) -pedantic -ffree-form
	$(G95)  $(FILES) -std=f2003 -pedantic -ffree-form

iverify:
	$(IFORT) -stand $(FS) -pedantic -ffree-form -warn stderrors $(FILES)

nverify:
	$(NFORT) $(FILES) -f$(FSTD) -free

pverify:
	$(PFORT) $(FILES) -Mfreeform

averify:
	$(AFORT) $(FILES) -f free

sverify:
	$(SFORT) $(FILES) -free

fverifym:
	$(FLANG)  $(M90FILES) -ffree-form
	cd $(INTERFORMDIR) && $(FLANG) $(MFILES) -fno-fixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(FLANG) $(MFILES) -fno-free-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

gverifym:
	$(GFORT)  $(M90FILES) -pedantic  -ffree-form
	cd $(INTERFORMDIR) && $(GFORT) $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(GFORT) $(MFILES) -pedantic -ffree-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

9verifym:
	$(G95) $(M90FILES) -pedantic -ffree-form
	cd $(INTERFORMDIR) && $(G95) $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(G95) $(MFILES) -pedantic -ffree-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

iverifym:
	$(IFORT) -standd:$(FS) -pedantic -ffree-form $(M90FILES)
	cd $(INTERFORMDIR) && $(IFORT) $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(IFORT) $(MFILES) -pedantic -ffree-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

nverifym:
	$(NFORT) $(M90FILES) -f$(FSTD) -free
	cd $(INTERFORMDIR) && $(NFORT) $(MFILES) -free && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(NFORT) ../common/consts.F ../common/debug.F ../common/memory.F ../common/info.F ../common/infnan.F ../common/lina.F ../common/output.F ../common/pintrf.F preproc.f initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f $(INTERSECTIONMEXDIR)/fmxapi.F $(INTERSECTIONMEXDIR)/prob.F $(INTERSECTIONMEXDIR)/classical/fmxcl.F $(INTERSECTIONMEXDIR)/gethuge.F $(INTERSECTIONMEXDIR)/newuoa-interface.F $(INTERSECTIONMEXDIR)/classical/newuoa-interface.F -I$(NEUPDFO) -fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

pverifym:
	$(PFORT) $(M90FILES) -Mfreeform
	cd $(INTERFORMDIR) && $(PFORT) $(MFILES) -Mfixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(PFORT) $(MFILES) -Mfreeform && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

averifym:
	$(AFORT) $(M90FILES) -f free
	cd $(INTERFORMDIR) && $(AFORT)  $(MFILES) -f fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(AFORT)  $(MFILES) -f free && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

sverifym:
	$(SFORT) $(M90FILES) -free
	cd $(INTERFORMDIR) && $(SFORT) $(MFILES) -fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)
	cd $(INTERFORMDIR) && $(SFORT) $(MFILES) -free && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

# Please change the variables FC and FFLAGS to match your fortran compiler

FC = sunf90
FFLAGS = -g
COMMON = ../common

# The .F90 or .f90 files will be compiled in the enumeration order of
# that the .o files. The order matters, because the compilation of
# each .o file depends on the .o files before it.

example: example.f90 \
	consts.o info.o debug.o memory.o infnan.o lina.o output.o pintrf.o \
	preproc.o initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o
	$(FC) $(FFLAGS) -o example example.f90 *.o

test: example
	./example

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

preproc.o: preproc.f90
	$(FC) $(FFLAGS) -c preproc.f90
initialize.o: initialize.f90
	$(FC) $(FFLAGS) -c initialize.f90
trustregion.o: trustregion.f90
	$(FC) $(FFLAGS) -c trustregion.f90
vlagbeta.o: vlagbeta.f90
	$(FC) $(FFLAGS) -c vlagbeta.f90
geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c geometry.f90
shiftbase.o: shiftbase.f90
	$(FC) $(FFLAGS) -c shiftbase.f90
update.o: update.f90
	$(FC) $(FFLAGS) -c update.f90
newuob.o: newuob.f90
	$(FC) $(FFLAGS) -c newuob.f90
newuoa.o: newuoa.f90
	$(FC) $(FFLAGS) -c newuoa.f90

clean:
	rm -f *.o *.mod a.out example *__genmod.f90 *.dbg
	cd $(COMMON) && rm -f *.o *.mod a.out example *__genmod.f90 *.dbg && cd $(CURRENTDIR)
