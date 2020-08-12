#GFORT=gfortran -g -Wall -Wextra -Wconversion -Warray-temporaries -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
GFORT=gfortran -g -Wall -Wextra -Wconversion -Wuninitialized -Wmaybe-uninitialized -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow 
IFORT=ifort -g -warn all -check all  -fpe0 -debug extended -fimplicit-none #-names as_is
G95=g95 -g -Wall -Wextra -fbounds-check -fimplicit-none 
NFORT=nagfor -g -Warn=allocation -Warn=constant_coindexing -Warn=reallocation -Warn=subnormal -colour=error:red,warn:magenta,info:yellow
PFORT=pgfortran  -g 
ABSOFT=af95  -g -m1 -en
SFORT=sunf95  -g -w3 -u -ansi -U

CURRENTDIR=$(PWD)
INTERFORMDIR=../intersection_form/newuoa/
MEXDIR= $(PWD)/../../matlab/mex_gateways/
INTERSECTIONMEXDIR= $(PWD)/../../matlab/mex_gateways/intersection_form/

COMFILES=../common/consts.F90 ../common/debug.F90 ../common/memory.F90 ../common/info.F90 ../common/infnan.F90 ../common/lina.F90 ../common/output.F90 ../common/pintrf.F90 
FILES= $(COMFILES) preproc.f90 initialize.f90 trustregion.f90 vlagbeta.f90 geometry.f90 update.f90 shiftbase.f90 newuob.f90 newuoa.f90 example.f90
 
M90FILES= $(FILES) $(MEXDIR)/fmxapi.F90 $(MEXDIR)/prob.F90 $(MEXDIR)/classical/fmxcl.F90 $(MEXDIR)/gethuge.F90 $(MEXDIR)/newuoa-interface.F90 $(MEXDIR)/classical/newuoa-interface.F90 -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/

MFILES= ../common/consts.F ../common/debug.F ../common/memory.F ../common/info.F ../common/infnan.F ../common/lina.F ../common/output.F ../common/pintrf.F preproc.f initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f example.f $(INTERSECTIONMEXDIR)/fmxapi.F $(INTERSECTIONMEXDIR)/prob.F $(INTERSECTIONMEXDIR)/classical/fmxcl.F $(INTERSECTIONMEXDIR)/gethuge.F $(INTERSECTIONMEXDIR)/newuoa-interface.F $(INTERSECTIONMEXDIR)/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/


FSTD=2003
FS=f03

verify: gverify iverify 9verify nverify pverify averify sverify
verifym: gverifym iverifym 9verifym nverifym pverifym averifym sverifym

gverify:
	$(GFORT)  $(FILES) -std=f$(FSTD) -pedantic -ffree-form

9verify:
	$(G95)  $(FILES) -std=f$(FSTD) -pedantic -ffree-form

iverify:
	$(IFORT) -stand:$(FS) -pedantic -ffree-form $(FILES) 

nverify:
	$(NFORT) $(FILES) -f$(FSTD) -free

pverify:
	$(PFORT) $(FILES) -Mfreeform

averify:
	$(ABSOFT) $(FILES) -f free

sverify:
	$(SFORT) $(FILES) -free


gverifym:
	$(GFORT) -c  $(M90FILES) -pedantic  -ffree-form
	cd $(INTERFORMDIR) && $(GFORT) -c $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

9verifym:
	$(G95) -c $(M90FILES) -pedantic -ffree-form
	cd $(INTERFORMDIR) && $(G95) -c $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

iverifym:
	$(IFORT) -c -standd:$(FS) -pedantic -ffree-form $(M90FILES) 
	cd $(INTERFORMDIR) && $(IFORT) -c $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

nverifym:
	$(NFORT) -c $(M90FILES) -f$(FSTD) -free
	cd $(INTERFORMDIR) && $(NFORT) -c $(MFILES) -fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

pverifym:
	$(PFORT) -c $(M90FILES) -Mfreeform
	cd $(INTERFORMDIR) && $(PFORT) -c $(MFILES) -Mfixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

averifym:
	$(ABSOFT) -c $(M90FILES) -f free
	cd $(INTERFORMDIR) && $(ABSOFT) -c  $(MFILES) -f fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

sverifym:
	$(SFORT) -c $(M90FILES) -free
	cd $(INTERFORMDIR) && $(SFORT) -c $(MFILES) -fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

clean:
	rm -f *.o *.mod a.out *__genmod.f90
