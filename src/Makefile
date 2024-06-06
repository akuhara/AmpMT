#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
FC     = mpif90
FFLAGS = -pg -Wall -pedantic -fbounds-check -O0 -std=f2008 \
         -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	 -fconvert=big-endian -g -fcheck=all

# Example 2. GNU fortran for fast computation
#FC     = mpif90
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check \
         -fconvert=big-endian

# Example 3. Intel compiler for fast computation
#FC      = mpiifort
#FFLAGS  = -assume byterecl -lmpi -xAVX -O3 -convert big_endian
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example A. General case
FFTW   = -I/usr/include -lfftw3 -L/usr/lib64

# Example B. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------
COV =
BINDIR = ../bin
CMT_BY_AMP = $(BINDIR)/cmt_by_amp

OBJS_CMT_BY_AMP = cmt_by_amp.o cls_radiation.o cls_observation.o cls_moment.o cls_model.o mod_random.o mod_mpi.o cls_parallel.o cls_mcmc.o cls_forward.o cls_output.o

TARGET = $(CMT_BY_AMP)
.PHONY: all
all: $(TARGET)

cmt_by_amp.o: cls_radiation.mod cls_observation.mod cls_moment.mod \
	      cls_model.mod mod_mpi.mod cls_mcmc.mod cls_parallel.mod cls_forward.mod cls_output.mod
cls_model.o: mod_random.mod
cls_prallel.o: mod_mpi.mod cls_mcmc.mod mod_random.mod
cls_mcmc.o: mod_random.mod cls_model.mod
cls_forward.o: cls_model.mod cls_moment.mod cls_observation.mod cls_radiation.mod
cls_output.o: cls_radiation.mod cls_moment.mod

$(CMT_BY_AMP): $(OBJS_CMT_BY_AMP)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $^ $(FFTW) -o $@


.PHONY: clean
clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) $(COV) $(FFTW) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
