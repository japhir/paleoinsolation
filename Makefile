MAKEFLAGS += --no-builtin-rules --no-builtin-variables

.PHONY: all clean
all: solution runsnvec insolation

solution: ems-plan3.dat
runsnvec: out.dat
orbpars: orb.dat
insolation: ins.dat


ems-plan3.dat:
	curl 'https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat' > "$@"

snvec_clone:
	git clone 'https://github.com/rezeebe/snvec'
	touch snvec_clone  # to track when this was copied in

snvec/snvec-3.7.5.c: snvec_clone snvec.patch
	patch "$@" < snvec.patch
	
snvec/snvec.x: snvec/snvec-3.7.5.c
	gcc -o snvec/snvec.x snvec/snvec-3.7.5.c -lm
	
out.dat: snvec/snvec.x
	./snvec/snvec.x -1e5 1 1 . ems-plan3.dat

ins.dat: out.dat paleoinsolation.f90.o
	./paleoinsolation.f90.o

# one program file
paleoinsolation.f90.o: paleoinsolation.f90 kind.f90.o data.f90.o interp.f90.o orb.f90.o insolation.f90.o shr_kind_mod.f90.o shr_const_mod.f90.o shr_log_mod.f90.o shr_orb_mod.f90.o 
	gfortran -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# my kind function defines a dp type according to best practices
kind.f90.o: kind.f90
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# my own functions
# readdata and writedata
data.f90.o: data.f90 kind.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# the locate function
interp.f90.o: interp.f90 kind.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# calculate insolation from orbital parameters
insolation.f90.o: insolation.f90 kind.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# functions that could be embedded into other packages:

# drop-in replacement for shr_orb_params in CDEPS
# https://github.com/ESCOMP/CDEPS/blob/main/share/shr_orb_mod.F90#L236
shr_orb_mod.f90.o: shr_orb_mod.F90 shr_kind_mod.f90.o shr_const_mod.f90.o shr_log_mod.f90.o data.f90.o interp.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# module dependencies for CDEPS
shr_kind_mod.f90.o: shr_kind_mod.F90
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# for some reason this is not compiling a .mod file...
#shr_infnan_mod.f90.o: shr_infnan_mod.F90.in
#	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

#shr_strconvert_mod.f90.o: shr_strconvert_mod.F90 shr_infnan_mod.f90.o
#	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

shr_log_mod.f90.o: shr_log_mod.F90 shr_kind_mod.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^
 
shr_const_mod.f90.o: shr_const_mod.F90 shr_kind_mod.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^


# drop-in replacement for orbpar in paleoToolkit
# https://github.com/CESM-Development/paleoToolkit/blob/master/PaleoCalAdjust/f90/modules/GISS_orbpar_subs.f90
orb.f90.o: orb.f90 kind.f90.o data.f90.o interp.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

clean:
	#-rm 'ems-plan3.dat'
	#-rm -rf 'snvec'
	#-rm 'snvec_clone'
	#-rm 'out.bin'
	#-rm 'out.dat'
	-rm 'ins.dat'
	-rm *.f90.o
	-rm *.mod
