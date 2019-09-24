#Makefile

#compiler
cf=gfortran

#compiler flags
cflags=-O -fbacktrace -fbounds-check #-fopenmp 

#program output name
pg=kepler_tst

#modules
modk=kind_m.o
mod0=const_m.o
mod1=trafo3_m.o
mod2=orbdyn_m.o
mod3=kepler_m.o
mod4=kepler_tst.o


obj= $(modk) $(mod0) $(mod1) $(mod2) $(mod3) $(mod4) # $(mod5) $(mod6) $(mod7) $(mod8) $(mod9) $(mod10)
targ:$(obj)
	$(cf) ${cflags} -o $(pg) $(obj)
	
%.o: %.f90
	$(cf) ${cflags} -c $<

%.o: %.f
	$(cf) ${cflags} -c $<
    
# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

