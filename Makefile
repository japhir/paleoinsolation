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

# directories
SRC_DIR := src
OBJ_DIR := build
MOD_DIR := $(OBJ_DIR)/mod
OUT_DIR := fout

# compiler flags
# FCFLAGS := -J$(MOD_DIR) -I$(MOD_DIR) -Isrc -Wall -Wextra -O2
FCFLAGS := -J$(MOD_DIR) -I$(MOD_DIR) -Isrc -Wall -Wextra -O2

# list of all source files
SRCS := $(wildcard $(SRC_DIR)/*.f90) $(wildcard $(SRC_DIR)/*.F90)

# source files
OBJS := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(SRCS:.f90=.o))
OBJS := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(OBJS:.F90=.o))

# main test/executable
MAIN := $(SRC_DIR)/paleoinsolation.f90
LIB := lib$(NAME).a
TEST_EXE := $(OUT_DIR)/paleoinsolation.exe

# declare all public targets
.PHONY: all solution buildsnvec runsnvec fortran insolation clean cleanall
all: insolation

# ------------------------------------------------------------
# Automatic Fortran deps
# ------------------------------------------------------------
DEPS := $(OBJ_DIR)/deps.mk

$(DEPS): $(SRCS)
	@mkdir -p $(OBJ_DIR) $(MOD_DIR)
	@echo "Generating Fortran dependencies → $(DEPS)"
	@echo "# Auto-generated Fortran module dependencies" > $(DEPS).tmp; \
	for f in $(SRCS); do \
	  srcbase=$$(basename $$f); \
	  obj="$(OBJ_DIR)/$${srcbase%.*}.o"; \
	  \
	  used_mods=$$( \
	    awk 'BEGIN{IGNORECASE=1} \
	         /^[[:space:]]*use[[:space:]]+/ { \
	           mod=$$2; gsub(/,.*|[()]/,"",mod); \
	           tolower(mod); \
	           if (mod !~ /^iso_/ && mod != "omp_lib" && mod != "ieee_arithmetic") print mod \
	         }' $$f | sort -u \
	  ); \
	  \
	  for mod in $$used_mods; do \
	    echo "$$obj: $(MOD_DIR)/$$mod.mod" >> $(DEPS).tmp; \
	  done; \
	  \
	  defined_mods=$$( \
	    awk 'BEGIN{IGNORECASE=1} \
	         /^[[:space:]]*module[[:space:]]+/ && $$2 != "procedure" { \
	           mod=$$2; gsub(/,.*/,"",mod); tolower(mod); print mod \
	         }' $$f \
	  ); \
	  for mod in $$defined_mods; do \
	    echo "$(MOD_DIR)/$$mod.mod: $$obj" >> $(DEPS).tmp; \
	  done; \
	done; \
	mv $(DEPS).tmp $(DEPS)

-include $(DEPS)

# 
$(LIB): $(DEPS) $(filter-out $(OBJ_DIR)/paleoinsolation.o,$(OBJS))
	@echo "Creating $@"
	$(AR) $@ $^

# create object files from Fortran source
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	@echo "Compiling $<"
	$(FC) $(FCFLAGS) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_DIR)
	@echo "Compiling $<"
	$(FC) $(FCFLAGS) -c -o $@ $<

# link and archive
$(TEST_EXE): $(LIB)
	@mkdir -p $(OUT_DIR)
	@echo "Linking $@"
	$(LD) -o $@ $(MAIN) $(FCFLAGS) $(LIB)

# compile the fortran routines
fortran: $(TEST_EXE)

# broad overview targets
solution: dat/ZB18a-plan3.dat dat/ZB20a-plan3.dat
buildsnvec: snvec/snvec.x
runsnvec: dat/PT-ZB18a_1-1.dat dat/PT-ZB20a_1-1.dat
insolation: fortran out/ZB18a_insolation.dat

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
	./fout/paleoinsolation.exe


# cleanup, filter to avoid removing source code by accident
clean:
	$(RM) -r $(OBJ_DIR) $(OUT_DIR) $(LIB)
	$(RM) $(wildcard out/*.dat)

cleanall: clean
	$(RM) $(wildcard dat/*.dat)
	$(RM) $(wildcard dat/*.bin)
	$(RM) -r 'snvec'
	$(RM) 'clonesnvec'
