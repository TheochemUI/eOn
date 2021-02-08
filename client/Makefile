TARGET_NAME=eonclient

#------------------------------------
#C compiler settings
CC=gcc
CXX=g++
LD=ld
AR=ar cru
RANLIB=ranlib
DEPFLAGS=-MMD

#if building omp lammps
#LDFLAGS += -fopenmp
#if building lammps with meam
#LDFLAGS += potentials/LAMMPS/libmeam.a
#LDFLAGS += potentials/LAMMPS/libreax.a
LDFLAGS += -lcurl -lpython3
CXXFLAGS += -g -Wall -Wno-unused-local-typedefs -Wno-unused-function -Wfatal-errors -lpython3 $(OPT) #-I flags depend on location of Python.h library
CXXFLAGS +=-w
DEBUG=yes

#CXXFLAGS += -Wall -Wno-unused-local-typedefs -Wno-unused-function -Wfatal-errors -lpython3 -I/export/tarpit/Python-3.7.4/ -I/export/tarpit/Python-3.7.4/Include/ $(OPT) #-I flags depend on location of Python.h library


OPT=-O0
ifndef DEBUG
    NDEBUG=1
    OPT=-O2
endif

#-----------------------------------
#MPI settings
ifdef EONMPIBGP
    EONMPI=1
endif

ifdef EONMPI
#    CC=mpicc
    CXX=mpicxx
    CXXFLAGS += -DEONMPI
    TARGET_NAME=eonclientmpi
endif

#------------------------------------
#Ensures all code is statically linked on a Linux machine
UNAME := $(shell uname)
MACHINE := $(shell uname -m)
ifeq ($(UNAME), Linux)
  LINUX=1
  ifeq ($(MACHINE), x86_64)
    LINUX64=1
  else
    LINUX32=1
  endif
  #MPI often doesn't mix well with static linking...
  ifndef EONMPI
	  #LDFLAGS += -static
  endif
endif
ifeq ($(UNAME), Darwin)
  OSX=1
endif
ifeq ($(UNAME), CYGWIN_NT-5.1)
  WIN32=1
endif

ifdef LINUX
    CXXFLAGS += -DLINUX
endif
ifdef OSX
    CXXFLAGS += -DOSX
endif
ifdef WIN32
    CXXFLAGS += -DWIN32
endif

#------------------------------------
#Fortran compiler settings
ifdef NO_FORTRAN
    CXXFLAGS += -DNO_FORTRAN
else
    # gfortran
    FC = $(CROSS)gfortran
    FFLAGS += $(OPT) -ffpe-trap=invalid,zero,overflow
    FAR = $(CROSS)ar cru
#    LDFLAGS += -L/usr/local/opt/lib/gcc/7 -lgfortran
    LDFLAGS += -L/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9 -lgfortran
    # ifort
    #FC = ifort
    #FFLAGS += $(OPT) -fpe0
    # Linux linking
    #FAR = ar cru
    #LDFLAGS += -limf -lsvml -lifcore -ldl
    # OSX linking
    #FAR = ifort -staticlib -o
endif

#------------------------------------
# define NDEBUG to remove debugging lines
ifdef NDEBUG
    DEBUG=No
    CXXFLAGS += -DNDEBUG
else
    CXXFLAGS += -g
    DEBUG=Yes
endif

#LDFLAGS += -lquadmath
# include Rules.mk
include buildRules.mk
