MAKEFLAGS += --no-builtin-rules --no-builtin-variables

.PHONY: all clean
all: solution runsnvec insolation

solution: ems-plan3.dat
runsnvec: out.dat
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
paleoinsolation.f90.o: paleoinsolation.f90 kind.f90.o data.f90.o insolation.f90.o 
	gfortran -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# several module files
kind.f90.o: kind.f90
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# this also creates insol.mod
# TODO: make this a result
insolation.f90.o: insolation.f90 kind.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

# this also creates data.mod
# TODO: make this a result
data.f90.o: data.f90 kind.f90.o
	gfortran -c -std=f2008 -ffree-form -g -fcheck=bounds -o $@ $^

clean:
	-rm 'ems-plan3.dat'
	-rm -rf 'snvec'
	-rm 'snvec_clone'
	-rm 'out.bin'
	-rm 'out.dat'
	-rm 'ins.dat'
	-rm *.f90.o
	-rm *.mod
