#GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Warray-temporaries -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Wuninitialized -Wmaybe-uninitialized -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow 
IFORT=ifort -c -g -warn all -check all  -fpe0 -debug extended -fimplicit-none
G95=g95 -c -g -Wall -Wextra -fbounds-check -fimplicit-none 
NFORT=nagfor -c -g
PFORT=pgfortran -c -g
ABSOFT=af95 -c -g -m1 -en

CURRENTDIR=$(PWD)
INTERFORMDIR=../intersection_form/newuoa
MEXDIR= $(PWD)/../../matlab/mex_gateways/

FILES= ../consts.F ../debug.F ../info.F ../infnan.F ../lina.F ../output.F ../prob.F initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f90 newuoa.f 
 
M90FILES= ../consts.F ../debug.F ../info.F ../infnan.F ../lina.F ../output.F ../prob.F initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f90 newuoa.f $(MEXDIR)/mexapi.F $(MEXDIR)/gethuge.F $(MEXDIR)/newuoa-interface.F $(MEXDIR)/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/

MFILES= ../consts.F ../debug.F ../info.F ../infnan.F ../lina.F ../output.F ../prob.F initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f $(MEXDIR)/mexapi.F $(MEXDIR)/gethuge.F $(MEXDIR)/newuoa-interface.F $(MEXDIR)/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/


FSTD=2003

gverify:
	$(GFORT)  $(FILES) -std=f$(FSTD) -pedantic -ffree-form

9verify:
	$(G95)  $(FILES) -std=f$(FSTD) -pedantic -ffree-form

iverify:
	$(IFORT) $(FILES) $(MFILES) -std=f$(FSTD) -pedantic

nverify:
	$(NFORT) $(FILES) $(MFILES) -f$(FSTD) 

pverify:
	$(PFORT) $(FILES) $(MFILES)

averify:
	$(ABSOFT) $(FILES) $(MFILES)


gverifym:
	$(GFORT)  $(M90FILES) -pedantic 
	cd $(INTERFORMDIR) && $(GFORT)  $(MFILES) -pedantic -ffixed-form && cd $(CURRENTDIR)

9verifym:
	$(G95)  $(M90FILES) -pedantic 
	cd $(INTERFORMDIR) && $(G95)  $(MFILES) -pedantic -ffixed-form && cd $(CURRENTDIR)

clean:
	rm -f *.o *.mod a.out *__genmod.f90
