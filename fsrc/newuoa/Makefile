#GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Warray-temporaries -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Wuninitialized -Wmaybe-uninitialized -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow 
IFORT=ifort -c -g -warn all -check all  -fpe0 -debug extended -fimplicit-none
G95=g95 -c -g -Wall -Wextra -fbounds-check -fimplicit-none 
NFORT=nagfor -c -g
PFORT=pgfortran -c -g

FILES=../consts.F ../debug.F ../info.F ../infnan.F ../lina.F ../output.F initialize.f trustregion.f vlagbeta.f geometry.f update.f shiftbase.f newuob.f newuoa.f 
MFILES= ../consts.F ../debug.F ../infnan.F ../info.F ../lina.F ../../matlab/mex_gateways/mexapi.F ../../matlab/mex_gateways/gethuge.F ../../matlab/mex_gateways/newuoa-interface.F ../../matlab/mex_gateways/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/



FSTD=2003

gverify:
	$(GFORT)  $(FILES) -std=f$(FSTD) -pedantic
#	$(GFORT)  $(FILES) -pedantic

9verify:
	$(G95)  $(FILES) -std=f$(FSTD) -pedantic

iverify:
	$(IFORT) $(FILES) -std=f$(FSTD) -pedantic

nverify:
	$(NFORT) $(FILES) -f$(FSTD) 

pverify:
	$(PFORT) $(FILES)


gverifym:
	$(GFORT)  $(FILES) $(MFILES) 

9verifym:
	$(G95)  $(FILES) $(MFILES) -pedantic

iverifym:
	$(IFORT)  $(FILES) $(MFILES) -std=f$(FSTD) -pedantic

nverifym:
	$(NFORT) $(FILES) $(MFILES) -f$(FSTD) 

pverifym:
	$(PFORT) $(FILES) $(MFILES)


gverifym:

clean:
	rm -f *.o *.mod a.out *__genmod.f90
