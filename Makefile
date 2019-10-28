# Makefile for F90 version of MINRESQLP for symmetric systems.
# To run this file:
#    make
# Or to remove object files and executable:
#    make clean
#
# Authors:
#     Sou-Cheng Choi <sctchoi@uchicago.edu>
#     Computation Institute (CI)
#     University of Chicago
#     Chicago, IL 60637, USA
#
#     Michael Saunders <saunders@stanford.edu>
#     Systems Optimization Laboratory (SOL)
#     Stanford University
#     Stanford, CA 94305-4026, USA
#
# 26 Aug 2012: First release version.
# 15 Oct 2007: Don't use compiler option -r8 to get double precision
#              because it is nonstandard.
#              Use minresqlpDataModule to define real(kind=dp).
# 11 Oct 2007: First version for compiling minresqlpTestProgram.f90
#              and associated modules.
#              All files listed are .f90 source code.
#              C is not currently used.

# Set exactly one of these to yes
# (Extend the list if necessary)

USEg95        = no
USEgfortran   = yes
USEgeneric90  = no
USEgeneric    = no
USEnag        = no

ifeq ($(USEg95),yes)
  FC      =  g95
  FFLAGS1 = -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full
  FFLAGS2 = -g -O
endif

ifeq ($(USEgfortran),yes)
  FC      =  gfortran
# Open MPI 1.7+ uses mpifort
  MPIFC   =  mpifort
  FFLAGS1 = -g -O0 -pedantic -Wall -W -fbounds-check
  FFLAGS2 = -g -O3
#
#  On AMD Epyc, add flags suggested by
#  http://www.prace-ri.eu/best-practice-guide-amd-epyc
#
ifeq ($(shell lscpu|grep -m 1 -c -i epyc),1)
FFLAGS += -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer
else
# Assume non-EPYC is Intel and use MKL:
  FFLAGS = -I/opt/intel/mkl/include -cpp
  LFLAGS = -L/opt/intel/mkl/lib/intel64/ -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
endif
endif

ifeq ($(USEgeneric),yes)
  FC      =  f95
  FFLAGS1 = -g -O0 -cpp
  FFLAGS2 = -g -O -cpp
endif

ifeq ($(USEgeneric90),yes)
  FC      =  f90
  FFLAGS1 = -g -O0
  FFLAGS2 = -g -O
endif

ifeq ($(USEnag),yes)       
  FC      =  nagfor
  FFLAGS1 = -C=all -C=undefined -nan -gline -f2003  -g90 -u -kind=byte
  FFLAGS2 = -nan -gline -f2003  -g90 -u -kind=byte
endif

# Select one of these
#FFLAGS  += ${FFLAGS1}    # for development
FFLAGS  += ${FFLAGS2}    # for release


%.o: %.f90
	${FC} ${FFLAGS} -c -o $@ $<
%_mpi.o: %.f90
	${MPIFC} ${FFLAGS} -DUSE_MPI -c -o $@ $<

files = minresqlpDataModule.o minresqlpBlasModule.o minresqlpModule.o mm_ioModule.o minresqlpReadMtxModule.o minresqlpTestModule.o minresqlpTestProgram.o

pfiles = minresqlpDataModule.o parallelBlasModule_mpi.o minresqlpBlasModule.o minresqlpModule_mpi.o mm_ioModule.o minresqlpReadMtxModule.o minresqlpTestModule_mpi.o minresqlpTestProgram_mpi.o

serial: minresqlptest
parallel: minresqlptest_mpi

minresqlptest: ${files}
	${FC} ${FFLAGS} ${LFLAGS} -o $@ ${files}

minresqlptest_mpi: ${pfiles}
	${MPIFC} ${FFLAGS} -lmpi ${LFLAGS} -o $@ ${pfiles}

clean:
	rm -f *.o *.mod minresqlptest minresqlptest_mpi
