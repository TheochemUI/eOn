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

#intel 12 linking
#LDFLAGS += -limf -lsvml -lifcore -ldl -lirc


CXXFLAGS += -Wall -Wno-unused-local-typedefs -Wno-unused-function -Wfatal-errors $(OPT)

OPT=-O0
ifndef DEBUG
    NDEBUG=1
    OPT=-O2
endif

#-----------------------------------
# Documentation Settings
DOXYCONFIG := ClientEON.dox
# BADPRACTICE: need to export variables so they are recognized
# by doxygen; there has to be some workaround, so this is a 
# temporary solution
export DOXYDIR := ./../doxygen
export VERSION := $(shell svnversion)
export BUILDDATE = $(shell date)
export BUILDHOST := $(shell hostname -d)
export BUILDUSER := $(shell whoami)

#-----------------------------------
#MPI settings
ifdef EONMPIBGP
    EONMPI=1
endif

ifdef EONMPI
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
    ifdef BOINC
        ifdef OSX
            LDFLAGS += /opt/local/lib/gcc45/libgfortran.a
        else
            LDFLAGS += -lgfortran
        endif
    else
        LDFLAGS += -lgfortran
    endif
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
# define BOINC to link against boinc libs
ifdef BOINC
    BOINC=Yes
    DEBUG=No
    OBJECTS += Compression.o
    CXXFLAGS += -DNDEBUG -DBOINC -Iboinc/include
    OBJECTS += Compression.o
    TARGET_NAME = eonclientboinc

    ifdef WIN32
        LDFLAGS += -static -Lboinc/lib_win -lboinc_api -lboinc 
        LDFLAGS += -Larchive/lib_win -larchive -lz 
        CXXFLAGS += -Iarchive/include -DLIBARCHIVE_STATIC
    endif

    ifdef OSX
        OSX_FLAGS=-mmacosx-version-min=10.4 -DMAC_OS_X_VERSION_MIN_REQUIRED=1040  -Wfatal-errors 
        CC=/usr/bin/gcc
        CXX=/usr/bin/g++
        CC += $(OSX_FLAGS)
        CXX += $(OSX_FLAGS)
        LDFLAGS += -Lboinc/lib_osx -lboinc_api -lboinc -lpthread -Iarchive/include archive/lib_osx/libarchive.a -lz $(OSX_FLAGS)
    endif

    ifdef LINUX
    UNAME := $(shell uname -m)
    ifeq ($(UNAME), x86_64)
        LDFLAGS += -static -Lboinc/lib_linux64 -lboinc_api -lboinc -lpthread -larchive -lz
    else
    LDFLAGS += -static -Lboinc/lib_linux32 -lboinc_api -lboinc -lpthread -larchive -lz
    endif
    endif
else
    BOINC=No
    # define NDEBUG to remove debugging lines
    ifdef NDEBUG
        DEBUG=No
        CXXFLAGS += -DNDEBUG
    else
        CXXFLAGS += -g
        DEBUG=Yes
    endif
endif
#LDFLAGS += -lquadmath
include Rules.mk
