#GFORT=gfortran -c -g -Wall -Wextra -Wconversion -Warray-temporaries -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
GFORT=gfortran -c -g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow 
IFORT=ifort -c -g -warn all -check all -fpe0 -debug extended -fimplicit-none
G95=g95 -c -g -Wall -Wextra -fbounds-check -fimplicit-none 

FILES=../consts.F ../warnerror.F ../infos.f ../infnan.F ../lina.F initialize.f trsubp.f geometry.f update.f calquad.f vlagbeta.f shiftbase.f newuob.f newuoa.f 
MFILES= ../consts.F ../warnerror.F ../infnan.F ../infos.f ../lina.F ../../matlab/mex_gateways/mex_mod.F ../../matlab/mex_gateways/gethuge.F ../../matlab/mex_gateways/newuoa-interface.F ../../matlab/mex_gateways/classical/newuoa-interface.F -I/home/zaikunzhang/local/MATLAB/R2018a/extern/include/



FSTD=2003

gverify:
	$(GFORT)  $(FILES) -std=f$(FSTD) -pedantic

9verify:
	$(G95)  $(FILES) -std=f$(FSTD) -pedantic

iverify:
	$(IFORT) $(FILES) -std=f$(FSTD) -pedantic

gverifym:
	$(GFORT)  $(FILES) $(MFILES) 

9verifym:
	$(G95)  $(FILES) $(MFILES) -pedantic

iverifym:
	$(IFORT)  $(FILES) $(MFILES) -std=f$(FSTD) -pedantic



clean:
	rm -f *.o *.mod a.out *__genmod.f90
