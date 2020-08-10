#GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Warray-temporaries -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Wuninitialized -Wmaybe-uninitialized -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow 
IFORT=ifort -c -g -warn all -check all  -fpe0 -debug extended -fimplicit-none
G95=g95 -c -g -Wall -Wextra -fbounds-check -fimplicit-none 
NFORT=nagfor -c -g
PFORT=pgfortran -c -g
ABSOFT=af95 -c -g -m1 -en

CURRENTDIR=$(PWD)
INTERFORMDIR=../intersection_form/newuoa/
MEXDIR= $(PWD)/../../matlab/mex_gateways/
INTERSECTIONMEXDIR= $(PWD)/../../matlab/mex_gateways/intersection_form/

COMFILES=../common/consts.F90 ../common/debug.F90 ../common/memory.F90 ../common/info.F90 ../common/infnan.F90 ../common/lina.F90 ../common/output.F90 ../common/prob.F90 
FILES= $(COMFILES) preproc.f90 initialize.f90 trustregion.f90 vlagbeta.f90 geometry.f90 update.f90 shiftbase.f90 newuob.f90 newuoa.f90 
 
M90FILES= $(FILES) $(MEXDIR)/fmxapi.F90 $(MEXDIR)/problem.F90 $(MEXDIR)/classical/fmxcl.F90 $(MEXDIR)/gethuge.F90 $(MEXDIR)/newuoa-interface.F90 $(MEXDIR)/classical/newuoa-interface.F90 -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/

MFILES= ../common/consts.F ../common/debug.F ../common/memory.F ../common/info.F ../common/infnan.F ../common/lina.F ../common/output.F ../common/prob.F preproc.f initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f $(INTERSECTIONMEXDIR)/fmxapi.F $(INTERSECTIONMEXDIR)/problem.F $(INTERSECTIONMEXDIR)/classical/fmxcl.F $(INTERSECTIONMEXDIR)/gethuge.F $(INTERSECTIONMEXDIR)/newuoa-interface.F $(INTERSECTIONMEXDIR)/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/


FSTD=2003
FS=f03

verify: gverify iverify 9verify nverify pverify averify
verifym: gverifym iverifym 9verifym nverifym pverifym averifym

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


gverifym:
	$(GFORT)  $(M90FILES) -pedantic  -ffree-form
	cd $(INTERFORMDIR) && $(GFORT)  $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

9verifym:
	$(G95)  $(M90FILES) -pedantic -ffree-form
	cd $(INTERFORMDIR) && $(G95)  $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

iverifym:
	$(IFORT) -standd:$(FS) -pedantic -ffree-form $(M90FILES) 
	cd $(INTERFORMDIR) && $(IFORT)  $(MFILES) -pedantic -ffixed-form && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

nverifym:
	$(NFORT) $(M90FILES) -f$(FSTD) -free
	cd $(INTERFORMDIR) && $(NFORT)  $(MFILES) -fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

pverifym:
	$(PFORT) $(M90FILES) -Mfreeform
	cd $(INTERFORMDIR) && $(PFORT)  $(MFILES) -Mfixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

averifym:
	$(ABSOFT) $(M90FILES) -f free
	cd $(INTERFORMDIR) && $(ABSOFT)  $(MFILES) -f fixed && rm -f *.o *.mod a.out *__genmod.f90 && cd $(CURRENTDIR)

clean:
	rm -f *.o *.mod a.out *__genmod.f90
