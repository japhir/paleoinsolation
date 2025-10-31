# fortran compile stuff
# disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# project name
NAME := paleoinsolation

# configuration settings
FC := gfortran
AR := ar rcs
LD := $(FC)
RM := rm -f

# list of all source files
SRCS := src/kind.f90 \
	src/data.f90 \
	src/interp.f90 \
	src/insolation.f90 \
	src/orb.f90 \
	src/shr_kind_mod.F90 \
	src/shr_log_mod.F90 \
	src/shr_strconvert_mod.F90 \
	src/shr_const_mod.F90 \
	src/shr_orb_mod.F90
TEST_SRCS := src/paleoinsolation.f90

# create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f90, %.exe, $(TEST_SRCS))

# declare all public targets
.PHONY: all solution buildsnvec runsnvec insolation clean cleanall
all: $(LIB) $(TEST_EXE) insolation

# compile the fortran routines

# create static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

# link the executable
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

# create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: %
	$(FC) -c -o $@ $<

# define all module interdependencies
kind.mod := src/kind.f90.o
data.mod := src/data.f90.o
interp.mod := src/interp.f90.o
orb.mod := src/orb.f90.o
insolation.mod := src/insolation.f90.o
paleoinsolation.mod := src/paleoinsolation.f90.o

shr_kind_mod.mod := src/shr_kind_mod.F90.o
shr_const_mod.mod := src/shr_const_mod.F90.o
shr_log_mod.mod := src/shr_log_mod.F90.o
shr_orb_mod.mod := src/shr_orb_mod.F90.o
shr_strconvert_mod.mod := src/shr_strconvert_mod.F90.o

src/data.f90.o: $(kind.mod)
src/interp.f90.o: $(kind.mod)
src/insolation.f90.o: $(kind.mod)

src/paleoinsolation.f90.o: $(kind.mod)
src/paleoinsolation.f90.o: $(data.mod)
src/paleoinsolation.f90.o: $(interp.mod)
src/paleoinsolation.f90.o: $(orb.mod)
src/paleoinsolation.f90.o: $(insolation.mod)

src/paleoinsolation.f90.o: $(shr_kind_mod.mod)
src/paleoinsolation.f90.o: $(shr_const_mod.mod)
src/paleoinsolation.f90.o: $(shr_log_mod.mod)
src/paleoinsolation.f90.o: $(shr_orb_mod.mod)

# drop-in replacement for shr_orb_params in CDEPS
# https://github.com/ESCOMP/CDEPS/blob/main/share/shr_orb_mod.F90#L236
src/shr_orb_mod.f90.o: $(data.mod)
src/shr_orb_mod.f90.o: $(interp.mod)
src/shr_orb_mod.f90.o: $(shr_kind.mod)
src/shr_orb_mod.f90.o: $(shr_const.mod)
src/shr_orb_mod.f90.o: $(shr_log.mod)
src/shr_log_mod.f90.o: $(shr_kind_mod.mod)
src/shr_const_mod.f90.o: $(shr_kind_mod.mod)

# drop-in replacement for orbpar in paleoToolkit
# https://github.com/CESM-Development/paleoToolkit/blob/master/PaleoCalAdjust/f90/modules/GISS_orbpar_subs.f90
src/orb.f90.o: $(kind.mod)
src/orb.f90.o: $(data.mod)
src/orb.f90.o: $(interp.mod)

# end of Fortran dependencies, remainder is snvec stuff and actually running the code

# broad overview targets
solution: dat/ZB18a-plan3.dat dat/ZB20a-plan3.dat
buildsnvec: snvec/snvec.x
runsnvec: dat/PT-ZB18a_1-1.dat dat/PT-ZB20a_1-1.dat
insolation: out/ZB18a_insolation.dat

### solution:
# download the orbital solution from the web
# this is the 100-Myr version
# we now provide the 300 Myr full files,
# for now embedded in this repo!
# ems-plan3.dat:
# 	curl 'https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/PrecTilt/OS/ZB18a/ems-plan3.dat' > "$@"

# extract the archive data files
ifeq ($(wildcard dat/ZB18a-plan3.dat),)
dat/ZB18a-plan3.dat:
	tar xvf 'dat/ZB18a-plan3.tar.gz' --directory 'dat'
else
dat/ZB18a-plan3.dat:
	echo "using local file $@"
endif

ifeq ($(wildcard dat/ZB20a-plan3.dat),)
dat/ZB20a-plan3.dat:
	tar xvf 'dat/ZB20a-plan3.tar.gz' --directory 'dat'
else
dat/ZB20a-plan3.dat:
	echo "using local file $@"
endif

# git clone the snvec c-program
clonesnvec:
	git clone 'https://github.com/rezeebe/snvec'
	touch clonesnvec # empty file to track when this was copied in

### buildsnvec:
# patch it to also save the lpx
snvec/snvec-3.7.5.c: clonesnvec snvec.patch
	patch "$@" < snvec.patch

# compile snvec
snvec/snvec.x: snvec/snvec-3.7.5.c
	gcc -std=c99 -o snvec/snvec.x snvec/snvec-3.7.5.c -lm

### runsnvec
# run snvec for ZB18a
ifeq ($(wildcard dat/PT-ZB18a_1-1.dat),)
dat/PT-ZB18a_1-1.dat: snvec/snvec.x dat/ZB18a-plan3.dat
	./snvec/snvec.x -3e5 1 1 'dat' 'ZB18a-plan3.dat'
	mv 'out.dat' 'dat/PT-ZB18a_1-1.dat'
	mv 'out.bin' 'dat/PT-ZB18a_1-1.bin'
else
dat/PT-ZB18a_1-1.dat:
	echo "using pre-built $@"
endif

# run snvec for ZB20a
ifeq ($(wildcard dat/PT-ZB20a_1-1.dat),)
dat/PT-ZB20a_1-1.dat: snvec/snvec.x dat/ZB20a-plan3.dat
	./snvec/snvec.x -3e5 1 1 'dat' 'ZB20a-plan3.dat'
	mv 'out.dat' 'dat/PT-ZB20a_1-1.dat'
	mv 'out.bin' 'dat/PT-ZB20a_1-1.bin'
else
dat/PT-ZB20a_1-1.dat:
	echo "using pre-built $@"
endif

ifeq ($(wildcard out),)
out:
	mkdir out
endif

### insolation
# run example fortran routine to calculate insolation
out/ZB18a_insolation.dat: out dat/PT-ZB18a_1-1.dat $(paleoinsolation.mod)
	./src/paleoinsolation.exe


# cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(TEST_EXE)
	$(RM) 'clonesnvec'
	$(RM) $(wildcard out/*.dat)
	$(RM) $(filter %.o, $(OBJS)) $(LIB) $(wildcard *.mod)

cleanall:
	$(RM) $(wildcard dat/*.dat)
	$(RM) $(wildcard dat/*.bin)
	$(RM) -r 'snvec'
